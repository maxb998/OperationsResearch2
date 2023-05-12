#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


static void setSEC(double *coeffs, int *indexes, CplexData *cpx, int *successors, int *subtoursMap, int subtourCount, int iterNum, Instance *inst, int nCols);

static void RepairHeuristic(int *successors, Solution *out, int *subtourMap, int subtourCount, double startTimeSec);

static inline void findBestSubtourMerge(Solution *sol, int mergeIndexes[2], int mergeIndexesSubtourIDs[2], int *invertOrientation, int subtourCount, size_t *subtourPosition);


Solution benders(Instance *inst, double tlim)
{
	size_t n = inst->nNodes;

	Solution sol1 = newSolution(inst), sol2 = newSolution(inst);
	Solution *best = &sol1, *repaired = &sol2;

    CplexData cpx = initCplexData(inst);

	struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	double currentTimeSec = startTimeSec;

	int ncols = CPXgetnumcols(cpx.env, cpx.lp);
	double *xstar = malloc(ncols * sizeof(double));

	int *successors = malloc(n * sizeof(int) * 2);
	int *subtoursMap = &successors[n];

	int *indexes = malloc(ncols * sizeof(int));

	while (currentTimeSec - startTimeSec < tlim)
	{
		static int iterNum = 0;

		// set time limit as remainig time from starting time
		clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    	currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
		CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, currentTimeSec + tlim - startTimeSec);

		if (CPXmipopt(cpx.env, cpx.lp))
			throwError(inst, NULL, "Benders: output of CPXmipopt != 0");

		if (CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1))
			throwError(inst, NULL, "Benders: output of CPXgetx != 0");

		int subtourCount = 0;
		
		cvtCPXtoSuccessors(xstar, ncols, n, successors, subtoursMap, &subtourCount);
		
		// generate a solution using Repair Heuristic and check if it is better than the previous solutions
		RepairHeuristic(successors, repaired, subtoursMap, subtourCount, startTimeSec);

		if (best->cost > repaired->cost)
		{
			Solution *temp;
			swapElems(best, repaired, temp);
		}

		LOG(LOG_LVL_LOG, "Number of subtours detected at iteration %d is %d", iterNum, subtourCount);

		if (subtourCount == 1) // means that there is only one subtour
			break;

		// add subtour elimination constraints
		double * coeffs = xstar; // reuse xStar instead of allocating new memory

		setSEC(coeffs, indexes, &cpx, successors, subtoursMap, subtourCount, iterNum, inst, ncols);
		
		iterNum++;
	}

	free(xstar);
	free(indexes);
	destroyCplexData(&cpx);

	destroySolution(repaired);
	free(successors);

    return *best;
}

static void setSEC(double *coeffs, int *indexes, CplexData *cpx, int *successors, int *subtoursMap, int subtourCount, int iterNum, Instance *inst, int nCols)
{
	size_t n = inst->nNodes;

	// set all coeffs to 1 at the beggining so we don't have to think about them again
	for (size_t i = 0; i < nCols; i++)
		coeffs[i] = 1.;

	static char sense = 'L';
	char *cname = malloc(20);
	static int izero = 0;

	for (int subtourID = 0; subtourID < subtourCount; subtourID++)
	{
		// get first node of next subtour
		int subtourStart = 0;
		for (;subtoursMap[subtourStart] < subtourID; subtourStart++);

		int nnz = 0;
		double rhs = -1;

		// follow successor and add all edges that connect each element of the subtour into the constraint
		int next = subtourStart;
		do
		{
			for (int i = successors[next]; i != subtourStart; i = successors[i])
			{
				indexes[nnz] = (int)xpos(next, i, n);
				nnz++;
			}
			rhs++;
			next = successors[next];
		} while (next != subtourStart);

		sprintf(cname, "SEC(%03d,%03d)", iterNum, subtourID);

		if (CPXaddrows(cpx->env, cpx->lp, 0, 1, nnz, &rhs, &sense, &izero, indexes, coeffs, NULL, &cname))
			throwError(inst, NULL, "benders: add SEC -> CPXaddrows(): error");
	}

	free(cname);
}

static void RepairHeuristic(int *successors, Solution *out, int *subtourMap, int subtourCount, double startTimeSec)
{
	Instance *inst = out->instance;
	size_t n = out->instance->nNodes;

	out->cost = 0.;

	// ##############################################################################
	// convert successor array to standard solution form adopted in Solution.indexpath

	// array containing the starting position of each subtour in the converted solution
	size_t *subtourPosition = malloc(subtourCount * sizeof(size_t));

	for (int subtourID = 0, pos = 0; subtourID < subtourCount; subtourID++)
	{
		// find first node of the subtour with id subtourID in successors

		subtourPosition[subtourID] = pos;

		size_t subtourFirstPos = 0;
		for (; subtourFirstPos < n; subtourFirstPos++)
			if (subtourMap[subtourFirstPos] == subtourID)
				break;
		
		// add elements to out.indexPath and .x and .y until we complete the subtour
		size_t i = subtourFirstPos;
		do
		{
			out->indexPath[pos] = i;
			out->X[pos] = inst->X[i];
			out->Y[pos] = inst->Y[i];
			pos++;
			i = successors[i];
		} while ((int)i != subtourFirstPos);
	}

	// ###############################################################################
	// Now to "repair" the solution
	while (subtourCount > 1)
	{
		// find the two subtours that merge with minimal cost increase (kind of like extra mileage heuristic but for subtours)

		// flag to check whether or not the direction of the loop might be inverted
		int invertOrientation;
		int mergeIndexes[2], mergeIndexesSubtourIDs[2];
		
		findBestSubtourMerge(out, mergeIndexes, mergeIndexesSubtourIDs, &invertOrientation, subtourCount, subtourPosition);
		
		size_t lastPos = n; // last position +1  of the of the subtour that will be merged and disappear (only used in the for condition)
		if (mergeIndexesSubtourIDs[1] < subtourCount-1)
			lastPos = subtourPosition[mergeIndexesSubtourIDs[1] + 1];

		size_t subtourSize = lastPos - subtourPosition[mergeIndexesSubtourIDs[1]];

		// create the "copy buffer" which contains the subtour rearranged in a way such that we just need to copy it into the correct position inside the other subtour

		double *xCopyBuff = malloc(subtourSize * (sizeof(double) * 2 + sizeof(int)));
		double *yCopyBuff = &xCopyBuff[subtourSize];
		int *indexCopyBuff = (int*)&yCopyBuff[subtourSize];
		
		size_t startPos = subtourPosition[mergeIndexesSubtourIDs[1]];
		for (size_t i = 0; i < subtourSize; i++)
		{
			size_t posToCopy;
			if (invertOrientation)
			{
				if (i < mergeIndexes[1] - startPos + 1)
					posToCopy = mergeIndexes[1] - i;
				else
					posToCopy = lastPos + (mergeIndexes[1] - startPos) - i;
			}
			else
			{
				if (i < lastPos - mergeIndexes[1] - 1)
					posToCopy = mergeIndexes[1] + 1 + i;
				else
					posToCopy = i - (lastPos - mergeIndexes[1] - 1) + startPos;
			}

			indexCopyBuff[i] = out->indexPath[posToCopy];
			xCopyBuff[i] = out->X[posToCopy];
			yCopyBuff[i] = out->Y[posToCopy];
		}

		// now to shift all data in order to make space inside subtour 0 to fit subtour 1 (0 and 1 are the subtour that will be merged)
		for (size_t i = subtourPosition[mergeIndexesSubtourIDs[1]]-1; i > mergeIndexes[0]; i--)
		{
			size_t shiftAmount = subtourSize;

			out->X[i + shiftAmount] = out->X[i];
			out->Y[i + shiftAmount] = out->Y[i];
			out->indexPath[i + shiftAmount] = out->indexPath[i];
		}

		// now to insert subtour 1 into subtour 0
		for (size_t i = 0; i < subtourSize; i++)
		{
			out->X[mergeIndexes[0] + 1 + i] = xCopyBuff[i];
			out->Y[mergeIndexes[0] + 1 + i] = yCopyBuff[i];
			out->indexPath[mergeIndexes[0] + 1 + i] = indexCopyBuff[i];
		}
		
		free(xCopyBuff);

		// fix auxiliary data structures
		subtourCount--;
		for (size_t i = mergeIndexesSubtourIDs[0] + 1; i < subtourCount; i++)
		{
			if (i < mergeIndexesSubtourIDs[1])
				subtourPosition[i] += subtourSize;
			else
				subtourPosition[i] = subtourPosition[i + 1];
		}
		
	}

	// set last position of the loop
	out->X[n] = out->X[0];
	out->Y[n] = out->Y[0];
	out->indexPath[n] = out->indexPath[0]; 

	out->cost = computeSolutionCostVectorized(out);

	if (inst->params.logLevel >= LOG_LVL_DEBUG)
		checkSolution(out);

	struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;

	out->execTime = currentTimeSec - startTimeSec;

	free(subtourPosition);
}

static inline void findBestSubtourMerge(Solution *sol, int mergeIndexes[2], int mergeIndexesSubtourIDs[2], int *invertOrientation, int subtourCount, size_t *subtourPosition)
{
	// THIS FUNCTION IS DEFINETIVELY VECTORIZABLE,
	// however I'm not going to do it since the extra code complexity isn't probably worth it since here the bottleneck is almost surely cplex's mipopt
	size_t n = sol->instance->nNodes;
	enum edgeWeightType ewt = sol->instance->params.edgeWeightType;
	int roundFlag = sol->instance->params.roundWeightsFlag;

	double min = INFINITY;

	for (size_t subtourID = 0; subtourID < subtourCount-1; subtourID++)
	{
		for (size_t subtourToCompare = subtourID+1; subtourToCompare < subtourCount; subtourToCompare++)
		{
			size_t last = subtourPosition[subtourToCompare+1];
			if (subtourToCompare == subtourCount - 1)
				last = n;

			for (size_t i = subtourPosition[subtourID]; i < subtourPosition[subtourID+1]; i++)
			{
				size_t secondI = i + 1;
				if (secondI == subtourPosition[subtourID+1])
					secondI = subtourPosition[subtourID];
				
				for (size_t j = subtourPosition[subtourToCompare]; j < last; j++ )
				{
					size_t secondJ = j + 1;
					if (secondJ == subtourPosition[subtourToCompare + 1])
						secondJ = subtourPosition[subtourToCompare];

					double cost = computeEdgeCost(sol->X[i], sol->Y[i], sol->X[secondJ], sol->Y[secondJ], ewt, roundFlag) +\
								   computeEdgeCost(sol->X[secondI], sol->Y[secondI], sol->X[j], sol->Y[j], ewt, roundFlag);
					
					if (cost < min)
					{
						min = cost;
						mergeIndexes[0] = i;
						mergeIndexes[1] = j;
						mergeIndexesSubtourIDs[0] = subtourID;
						mergeIndexesSubtourIDs[1] = subtourToCompare;
						*invertOrientation = 0;
					}

					cost = computeEdgeCost(sol->X[i], sol->Y[i], sol->X[j], sol->Y[j], ewt, roundFlag) +\
						   computeEdgeCost(sol->X[secondI], sol->Y[secondI], sol->X[secondJ], sol->Y[secondJ], ewt, roundFlag);

					if (cost < min)
					{
						min = cost;
						mergeIndexes[0] = i;
						mergeIndexes[1] = j;
						mergeIndexesSubtourIDs[0] = subtourID;
						mergeIndexesSubtourIDs[1] = subtourToCompare;
						*invertOrientation = 1;
					}
				}
			}
		}
	}
}

