#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <stdio.h>
#include <stdarg.h> // used for logger va_list



static inline void findBestSubtourMerge(SubtoursData *sub, int subtoursCount, Instance *inst, int mergeIndexes[2], int *invertOrientation);

static double computeSuccessorsSolCost(int *successors, Instance *inst);

CplexData initCplexData(Instance *inst)
{
	CplexData cpxData = { 0 };

	int errno = 0;
	cpxData.env = CPXopenCPLEX(&errno);
	if (errno)
		throwError(inst, NULL, "buildCPXModel: error at CPXopenCPLEX with code %d", errno);

	// screen output
	if (inst->params.logLevel >= LOG_LVL_EVERYTHING)
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
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
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
				cplexError(&cpxData, inst, NULL, "initCplexData: wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(cpxData.env, cpxData.lp)-1 != xpos(i,j,n) )
				cplexError(&cpxData, inst, NULL, "initCplexData: wrong position for x var.s");
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
			cplexError(&cpxData, inst, NULL, "initCplexData: CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	return cpxData;
}

void destroyCplexData(CplexData * cpxData)
{
	CPXfreeprob(cpxData->env, &cpxData->lp);
    CPXcloseCPLEX(&cpxData->env);
}

void cplexError(CplexData *cpxData, Instance *inst, Solution *sol, char *line, ...)
{
	printf("\033[1;31mERR \033[0m");

	va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);

    printf("\n");

    // free allocated memory
    if (inst) destroyInstance(inst);
    if (sol) destroySolution(sol);
	if (cpxData) destroyCplexData(cpxData);

    exit(EXIT_FAILURE);
}

int callbackError(char *line, ...)
{
	printf("\033[1;31mERR \033[0m");

	va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);

    printf("\n");

    exit(EXIT_FAILURE);
}

size_t xpos(size_t i, size_t j, size_t n)
{
    if (i == j)
        throwError(NULL, NULL, "xpos: i == j");
    if (i > j)
    {
        register size_t temp;
        swapElems(i, j, temp);
    }

    int pos = i * n + j - ((i + 1) * (i + 2)) / 2;
    return pos;
}

void cvtCPXtoSuccessors(double *xstar, int ncols, size_t nNodes, SubtoursData *subData)
{
	// reset arrays
	for (size_t i = 0; i < nNodes; i++)
		subData->successors[i] = -1;
	for (size_t i = 0; i < nNodes; i++)
		subData->subtoursMap[i] = -1;

	LOG(LOG_LVL_EVERYTHING, "###################################################");

	for (size_t i = 0; i < nNodes; i++)
	{
		size_t succ = i;
		for (size_t j = 0; j < nNodes; j++)
		{
			if ((succ != j) && (xstar[xpos(succ, j, nNodes)] > 0.5) && (subData->subtoursMap[j] == -1))
			{
				subData->successors[succ] = (int)j;
				subData->subtoursMap[succ] = subData->subtoursCount;
				LOG(LOG_LVL_EVERYTHING, "x(%3d,%3d) = 1   subtour n° %d\n", succ, j, subData->subtoursMap[succ] + 1);
				succ = j;
				j = 0;
			}
		}
		if (succ != i)
		{
			subData->successors[succ] = i;
			subData->subtoursMap[succ] = subData->subtoursCount;
			LOG(LOG_LVL_EVERYTHING, "x(%3d,%3d) = 1   subtour n° %d\n", succ, i, subData->subtoursMap[succ] + 1);
			(subData->subtoursCount)++;
		}
	}
}

void cvtSuccessorsToSolution(int *successors, Solution *sol)
{
	Instance *inst = sol->instance;
	size_t n = inst->nNodes;
	
	// start from 0
	sol->indexPath[0] = 0;
	sol->X[0] = inst->X[0];
	sol->Y[0] = inst->Y[0];

	for (size_t i = successors[0], pos = 1; i != 0; i = successors[i])
	{
		sol->indexPath[pos] = i;
		sol->X[pos] = inst->X[i];
		sol->Y[pos] = inst->Y[i];

		pos++;
	}

	sol->indexPath[n] = 0;
	sol->X[n] = inst->X[0];
	sol->Y[n] = inst->Y[0];
}

void setSEC(double *coeffs, int *indexes, CplexData *cpx, CPXCALLBACKCONTEXTptr context, SubtoursData *subData, int iterNum, Instance *inst, int nCols, int isBenders)
{
	size_t n = inst->nNodes;

	// set all coeffs to 1 at the beggining so we don't have to think about them again
	for (size_t i = 0; i < nCols; i++)
		coeffs[i] = 1.;

	static char sense = 'L';
	char *cname = malloc(20);
	static int izero = 0;

	for (int subtourID = 0; subtourID < subData->subtoursCount; subtourID++)
	{
		// get first node of next subtour
		int subtourStart = 0;
		for (;subData->subtoursMap[subtourStart] < subtourID; subtourStart++);

		int nnz = 0;
		double rhs = -1;

		// follow successor and add all edges that connect each element of the subtour into the constraint
		int next = subtourStart;
		do
		{
			for (int i = subData->successors[next]; i != subtourStart; i = subData->successors[i])
			{
				indexes[nnz] = (int)xpos(next, i, n);
				nnz++;
			}
			rhs++;
			next = subData->successors[next];
		} while (next != subtourStart);

		sprintf(cname, "SEC(%03d,%03d)", iterNum, subtourID);
		if(isBenders)
		{
			if (CPXaddrows(cpx->env, cpx->lp, 0, 1, nnz, &rhs, &sense, &izero, indexes, coeffs, NULL, &cname))
				cplexError(cpx, inst, NULL, "setSEC ->benders: CPXaddrows() error");
		}
		else
		{
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, indexes, coeffs))
				cplexError(cpx, inst, NULL, "setSEC ->lazyCallback: CPXaddlazyconstraints() error");
		}
	}

	free(cname);
}

double RepairHeuristicSuccessors(SubtoursData *sub, Instance *inst)
{
	size_t n = inst->nNodes;

	int subtoursCount = sub->subtoursCount;
	while (subtoursCount > 1)
	{
		int invert = 0;
		int mergeIndexes[2] = {0};

		findBestSubtourMerge(sub, subtoursCount, inst, mergeIndexes, &invert);

		if (invert == 0)
		{
			register int temp = sub->successors[mergeIndexes[0]];
			sub->successors[mergeIndexes[0]] = sub->successors[mergeIndexes[1]];
			sub->successors[mergeIndexes[1]] = temp;
		}
		else
		{
			int last = sub->successors[mergeIndexes[1]];
			int previous = sub->successors[mergeIndexes[1]];
			int current = sub->successors[previous];
			int next = sub->successors[current];
			do
			{
				sub->successors[current] = previous;
				previous = current;
				current = next;
				next = sub->successors[next];
			} while (current != last);
			
			sub->successors[last] = sub->successors[mergeIndexes[0]];
			sub->successors[mergeIndexes[0]] = mergeIndexes[1];
		}

		// now adjust subtoursMap
		int sub2Map = sub->subtoursMap[mergeIndexes[1]];
		for (size_t i = 0; i < n; i++)
		{
			if (sub->subtoursMap[i] == sub2Map)
				sub->subtoursMap[i] = sub->subtoursMap[mergeIndexes[0]];
			else if (sub->subtoursMap[i] > sub2Map)
				sub->subtoursMap[i]--;
		}

		subtoursCount--;
	}

	double cost = computeSuccessorsSolCost(sub->successors, inst);

	return cost;
}

static inline void findBestSubtourMerge(SubtoursData *sub, int subtoursCount, Instance *inst, int mergeIndexes[2], int *invertOrientation)
{
	float *X = inst->X, *Y =inst->Y;
	enum EdgeWeightType ewt = inst->params.edgeWeightType ;
	int roundFlag = inst->params.roundWeightsFlag;

	float min = INFINITY;

	for (size_t subtourID = 0; subtourID < subtoursCount-1; subtourID++)
	{
		size_t first1 = 0;
			while (sub->subtoursMap[first1] != subtourID)
				first1++;

		for (size_t subtourToCompare = subtourID+1; subtourToCompare < subtoursCount; subtourToCompare++)
		{
			size_t first2 = 0;
			while (sub->subtoursMap[first2] != subtourToCompare)
				first2++;
			
			int finish1 = 0;
			int i = first1;
			while ((!finish1) || (i != sub->successors[first1]))
			{
				if (sub->successors[i] == first1) finish1 = -1;
				
				// successor of i'th node
				int succI = sub->successors[i];

				int finish2 = 0;
				int j = first2;
				while ((!finish2) || (j != sub->successors[first2]))
				{
					if (sub->successors[j] == first2) finish2 = -1;

					// successor of j'th node
					int succJ = sub->successors[j];

					float cost = computeEdgeCost(X[i], Y[i], X[succJ], Y[succJ], ewt, roundFlag) + computeEdgeCost(X[succI], Y[succI], X[j], Y[j], ewt, roundFlag);

					if (cost < min)
					{
						min = cost;
						mergeIndexes[0] = i; mergeIndexes[1] = j;
						*invertOrientation = 0;
					}
					
					cost = computeEdgeCost(X[i], Y[i], X[j], Y[j], ewt, roundFlag) + computeEdgeCost(X[succI], Y[succI], X[succJ], Y[succJ], ewt, roundFlag);

					if (cost < min)
					{
						min = cost;
						mergeIndexes[0] = i; mergeIndexes[1] = j;
						*invertOrientation = 1;
					}

					j = succJ;
				}
				i = succI;
			}
		}
	}
}

static double computeSuccessorsSolCost(int *successors, Instance *inst)
{
	int n = (int)inst->nNodes;
	enum EdgeWeightType ewt = inst->params.edgeWeightType ;
	int roundFlag = inst->params.roundWeightsFlag;

	double cost = 0.0;

	int i = 0;
	int counter = 0;
	do
	{
		int succ = successors[i];
		cost += computeEdgeCost(inst->X[i], inst->Y[i], inst->X[succ], inst->Y[succ], ewt, roundFlag);
		i = succ;

		if (counter > n)
			throwError(inst, NULL, "computeSuccessorsSolCost: There are subtours inside the successor array even after repair heuristic");
		counter++;
	} while (i != 0);
	
	return cost;
}

int checkSuccessorSolution(Instance *inst, int *successors)
{
	size_t n = inst->nNodes;

	int *indexChecked = malloc(n * sizeof(int));
	indexChecked[0] = -1;
	for (size_t i = 1; i < n; i++)
		indexChecked[i] = 0;

	int i = 0, count = 0;
	do
	{
		if (indexChecked[i] == 1)
			return 1;
		indexChecked[i]++;
		i = successors[i];
		count++;
	} while (count <= n);
	
	return 0;
}

