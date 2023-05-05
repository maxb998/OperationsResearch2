#include "TspCplex.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

typedef struct
{
	CPXENVptr env;
	CPXLPptr lp;
	Instance *inst;
} CplexData;


static size_t xpos(size_t i, size_t j, size_t n);

static CplexData initCplexData(Instance *inst);

static void destroyCplexData(CplexData * cpxData);



static void RepairHeuristic(int *successors, Solution *out, int *subtourMap, int subtourCount);

static inline void findBestSubtourMerge(Solution *sol, int mergeIndexes[2], int mergeIndexesSubtourIDs[2], int *invertOrientation, int subtourCount, size_t *subtourPosition);


Solution blenders(Instance *inst, double tlim)
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

	int *successors = malloc(n * sizeof(int));
	int *subtoursMap = malloc(n * sizeof(int));

	int *indexes = malloc(ncols * sizeof(int));

	while (currentTimeSec - startTimeSec < tlim)
	{
		static int iterNum = 0;

		// set time limit as remainig time from starting time
		clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    	currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
		CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, currentTimeSec + tlim - startTimeSec);

		if (CPXmipopt(cpx.env, cpx.lp))
			throwError(inst, NULL, "Blenders: output of CPXmipopt != 0");

		if (CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1))
			throwError(inst, NULL, "Blenders: output of CPXgetx != 0");

		int subtourCount = 0;
		// reset arrays
		for (size_t i = 0; i < n; i++)
			successors[i] = -1;
		for (size_t i = 0; i < n; i++)
			subtoursMap[i] = -1;
		
		for (size_t i = 0; i < n; i++)
		{
			size_t succ = i;
			for (size_t j = 0; j < n; j++)
			{
				if ((succ != j) && (xstar[xpos(succ, j, n)] > 0.5) && (subtoursMap[j] == -1))
				{
					successors[succ] = (int)j;
					subtoursMap[succ] = subtourCount;
					LOG(LOG_LVL_EVERYTHING, "x(%3d,%3d) = 1   subtour n° %d\n", succ, j, subtoursMap[succ]);
					succ = j;
					j = 0;
				}
			}
			if (succ != i)
			{
				successors[succ] = i;
				subtoursMap[succ] = subtourCount;
				LOG(LOG_LVL_EVERYTHING, "x(%3d,%3d) = 1   subtour n° %d\n", succ, i, subtoursMap[succ]);
				subtourCount++;
			}
		}

		// generate a solution using Repair Heuristic and check if it is better than the previous solutions
		RepairHeuristic(successors, repaired, subtoursMap, subtourCount);
		if (best->cost > repaired->cost)
		{
			Solution *temp;
			swapElems(best, repaired, temp);
		}

		LOG(LOG_LVL_DEBUG, "Number of subtours detected at iteration %d is %d", iterNum, subtourCount);

		if (subtourCount == 1) // means that there is only one subtour
			break;

		// add subtour elimination constraints
		double * coeffs = xstar; // reuse xStar instead of allocating new memory

		for (int subtourID = 1; subtourID < subtourCount; subtourID++)
		{
			static char sense = 'L';
			char *cname = malloc(20);
			sprintf(cname, "SEC(%03d,%03d)", iterNum, subtourID);
			int nnz = 0;

			for (size_t i = 0; i < n; i++)
			{
				if (subtoursMap[i] == subtourID)
				{
					coeffs[nnz] = 1.;
					indexes[nnz] = (int)xpos(i, successors[i], n);
					nnz++;
				}
			}

			double rhs = nnz-1;
			static int izero = 0;

			if (CPXaddrows(cpx.env, cpx.lp, 0, 1, nnz, &rhs, &sense, &izero, indexes, coeffs, NULL, &cname))
				throwError(inst, NULL, "blenders: add SEC -> CPXaddrows(): error");

			free(cname);
		}
		
		iterNum++;
	}

	free(xstar);
	free(indexes);
	destroyCplexData(&cpx);

	destroySolution(repaired);
	free(subtoursMap);
	free(successors);

    return *best;
}

static void RepairHeuristic(int *successors, Solution *out, int *subtourMap, int subtourCount)
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

	out->cost = computeSolutionCostVectorizedDouble(out);

	if (inst->params.logLevel >= LOG_LVL_DEBUG)
		checkSolution(out);

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

static size_t xpos(size_t i, size_t j, size_t n)
{ 
	if ( i == j )
		throwError(NULL, NULL,"xpos: i == j");
	if ( i > j )
	{
		register size_t temp;
		swapElems(i,j,temp);
	}

	int pos = i * n + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}

static CplexData initCplexData(Instance *inst)
{
	CplexData cpxData;

	int errno = 0;
	cpxData.env = CPXopenCPLEX(&errno);
	if (errno)
		throwError(inst, NULL, "buildCPXModel: error at CPXopenCPLEX with code %d", errno);

	// screen output
	CPXsetintparam(cpxData.env, CPX_PARAM_SCRIND, CPX_ON);

	// random seed for cplex
	if (inst->params.randomSeed != -1)
		CPXsetintparam(cpxData.env, CPX_PARAM_RANDOMSEED, inst->params.randomSeed);
	else
		CPXsetintparam(cpxData.env, CPX_PARAM_RANDOMSEED, time(NULL));

	cpxData.lp = CPXcreateprob(cpxData.env, &errno, "TSP");
	if (errno)
		throwError(inst, NULL, "buildCPXModel: error at CPXcreateprob with code %d", errno);

    size_t n = inst->nNodes;
    enum edgeWeightType ewt = inst->params.edgeWeightType;
    int roundFlag = inst->params.roundWeightsFlag;

	char binary = CPX_BINARY; 

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

    // add binary var.s x(i,j) for i < j  
	
	for ( int i = 0; i < n; i++ )
	{
		for ( int j = i+1; j < n; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = computeEdgeCost(inst->X[i], inst->Y[i], inst->X[j], inst->Y[j], ewt, roundFlag); // cost == distance
			double ub = 1.0;
			if ( CPXnewcols(cpxData.env, cpxData.lp, 1, &obj, NULL, &ub, &binary, cname) ) 
				throwError(inst, NULL, "initCplexData: wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(cpxData.env, cpxData.lp)-1 != xpos(i,j,n) )
				throwError(inst, NULL, "initCplexData: wrong position for x var.s");
		}
	} 

    // add the degree constraints 

	int *index = (int *) calloc(n, sizeof(int));
	double *value = (double *) calloc(n, sizeof(double));

	for ( size_t h = 0; h < n; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%lu)", h+1);
		int nnz = 0;
		for ( size_t i = 0; i < n; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h,n);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(cpxData.env, cpxData.lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) 
			throwError(inst, NULL, "initCplexData: CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	return cpxData;
}

static void destroyCplexData(CplexData * cpxData)
{
	CPXfreeprob(cpxData->env, &cpxData->lp);
    CPXcloseCPLEX(&cpxData->env);
}
