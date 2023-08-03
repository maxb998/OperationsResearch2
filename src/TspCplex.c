#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <stdio.h>
#include <stdarg.h> // used for logger va_list


CplexData initCplexData(Instance *inst)
{
	CplexData cpxData = { .inst = inst };

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

    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

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
			{
				destroyCplexData(&cpxData);
				throwError(inst, NULL, "initCplexData: wrong CPXnewcols on x var.s");
			}
    		if ( CPXgetnumcols(cpxData.env, cpxData.lp)-1 != xpos(i,j,n) )
			{
				destroyCplexData(&cpxData);
				throwError(inst, NULL, "initCplexData: wrong position for x var.s");
			}
		}
	} 

    // add the degree constraints 

	int *index = (int *) calloc(n, sizeof(int));
	double *value = (double *) calloc(n, sizeof(double));

	for ( int h = 0; h < n; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);
		int nnz = 0;
		for ( int i = 0; i < n; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h,n);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(cpxData.env, cpxData.lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) )
		{
			destroyCplexData(&cpxData);
			throwError(inst, NULL, "initCplexData: CPXaddrows(): error 1");
		}
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

int xpos(int i, int j, int n)
{
    if (i == j)
        throwError(NULL, NULL, "xpos: i == j");
    if (i > j)
    {
        register int temp;
        swapElems(i, j, temp);
    }

    int pos = i * n + j - ((i + 1) * (i + 2)) / 2;
    return pos;
}

void cvtCPXtoSuccessors(double *xstar, int ncols, int nNodes, SubtoursData *subData)
{
	// reset arrays
	for (int i = 0; i < nNodes; i++)
		subData->successors[i] = -1;
	for (int i = 0; i < nNodes; i++)
		subData->subtoursMap[i] = -1;

	LOG(LOG_LVL_EVERYTHING, "###################################################");

	for (int i = 0; i < nNodes; i++)
	{
		int succ = i;
		for (int j = 0; j < nNodes; j++)
		{
			if ((succ != j) && (xstar[xpos(succ, j, nNodes)] > 0.5) && (subData->subtoursMap[j] == -1))
			{
				subData->successors[succ] = j;
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
	int n = inst->nNodes;

	if (sol->instance->params.logLevel >= LOG_LVL_DEBUG)
		checkSuccessorSolution(sol->instance, successors);
	
	// start from 0
	sol->indexPath[0] = 0;

	for (int i = successors[0], pos = 1; i != 0; i = successors[i], pos++)
		sol->indexPath[pos] = i;

	sol->indexPath[n] = 0;

	sol->cost = computeSolutionCost(sol);
}

void cvtSolutionToSuccessors(Solution *sol, int* successors)
{
	int n = sol->instance->nNodes;

	for (int i = 0; i < n; i++)
		successors[sol->indexPath[i]] = sol->indexPath[i+1];

	if ((sol->instance->params.logLevel >= LOG_LVL_EVERYTHING) && (!checkSuccessorSolution(sol->instance, successors) != 0))
		throwError(sol->instance, sol, "cvtSolutionToSuccessors: Converted solution is wrong");	
}

int setSEC(double *coeffs, int *indexes, CplexData *cpx, CPXCALLBACKCONTEXTptr context, SubtoursData *subData, int iterNum, Instance *inst, int nCols, bool isBenders)
{
	int retVal = 0;
	int n = inst->nNodes;

	// set all coeffs to 1 at the beggining so we don't have to think about them again
	for (int i = 0; i < nCols; i++)
		coeffs[i] = 1.;

	static char sense = 'L';
	char *cname = malloc(20);
	static int izero = 0;

	for (int subtourID = 0; (subtourID < subData->subtoursCount) && (retVal == 0); subtourID++)
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
				indexes[nnz] = xpos(next, i, n);
				nnz++;
			}
			rhs++;
			next = subData->successors[next];
		} while (next != subtourStart);

		sprintf(cname, "SEC(%03d,%03d)", iterNum, subtourID);

		if(isBenders)
			retVal = CPXaddrows(cpx->env, cpx->lp, 0, 1, nnz, &rhs, &sense, &izero, indexes, coeffs, NULL, &cname);
		else
			retVal = CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, indexes, coeffs);
	}

	free(cname);

	return retVal;
}

__uint128_t computeSuccessorsSolCost(int *successors, Instance *inst)
{
	if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (!checkSuccessorSolution(inst, successors) != 0))
		throwError(inst, NULL, "computeSuccessorsSolCost: successors array passed as input is not feasible");

	int n = inst->nNodes;
	enum EdgeWeightType ewt = inst->params.edgeWeightType ;
	bool roundFlag = inst->params.roundWeights;

	__uint128_t cost = 0LL;

	int i = 0;
	int counter = 0;
	do
	{
		int succ = successors[i];
		cost += cvtFloat2Cost(computeEdgeCost(inst->X[i], inst->Y[i], inst->X[succ], inst->Y[succ], ewt, roundFlag));
		i = succ;

		if (counter > n)
			throwError(inst, NULL, "computeSuccessorsSolCost: There are subtours inside the successor array even after repair heuristic");
		counter++;

	} while (i != 0);
	
	return cost;
}

bool checkSuccessorSolution(Instance *inst, int *successors)
{
	int n = inst->nNodes;

	int *indexChecked = malloc(n * sizeof(int));
	indexChecked[0] = -1;
	for (int i = 1; i < n; i++)
		indexChecked[i] = 0;

	bool retVal = true;

	int i = 0, count = 0;
	do
	{
		if (indexChecked[i] == 1)
		{ LOG(LOG_LVL_CRITICAL, "checkSuccessorSolution: Node %d appears multiple times in successors", i); retVal = false; break; }

		indexChecked[i]++;
		i = successors[i];
		count++;
	} while (count <= n);
	
	return retVal;
}

int WarmStart(CplexData *cpx, int *successors)
{
	Instance *inst = cpx->inst;
	int n = inst->nNodes;

	if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (!checkSuccessorSolution(inst, successors)))
	{
		LOG(LOG_LVL_ERROR, "WarStartSuccessors: successor solution is incorrect");
		return 1;
	}

	double *ones = malloc(n * sizeof(double));
	int *indexes = malloc(n * sizeof(int));
	for (int i = 0; i < n; i++)
		ones[i] = 1.0;
	for (int i = 0; i < n; i++)
		indexes[i] = xpos(i, successors[i], n);

	int izero = 0;
	char *mipstartName = "Warm Start, External";
	int effort = CPX_MIPSTART_NOCHECK;
	if (inst->params.logLevel >= LOG_LVL_DEBUG)
		effort = CPX_MIPSTART_CHECKFEAS;

	int retVal = CPXaddmipstarts(cpx->env, cpx->lp, 1, n, &izero, indexes, ones, &effort, &mipstartName);

	free(ones);
	free(indexes);

	return retVal;
}


SubtoursData initSubtoursData(int n)
{
	SubtoursData subData = {
		.subtoursCount = 0,
		.successors = malloc((n + AVX_VEC_SIZE) * sizeof(int)),
		.subtoursMap = malloc((n + AVX_VEC_SIZE) * sizeof(int))
	};

	if ((subData.successors == NULL) || (subData.subtoursMap == NULL))
		throwError(NULL, NULL, "initSubtoursData: Failed to allocate memory");

	return subData;
}


void destroySubtoursData(SubtoursData *subData)
{
	free(subData->successors);
	free(subData->subtoursMap);
	subData->successors = NULL;
	subData->subtoursMap = NULL;
}
