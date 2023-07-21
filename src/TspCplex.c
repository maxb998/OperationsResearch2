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

    size_t n = inst->nNodes;
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

void cvtSolutionToSuccessors(Solution *sol, int* successors)
{
	size_t n = sol->instance->nNodes;

	for (size_t i = 0; i < n; i++)
		successors[sol->indexPath[i]] = sol->indexPath[i+1];

	if ((sol->instance->params.logLevel >= LOG_LVL_EVERYTHING) && (!checkSuccessorSolution(sol->instance, successors) != 0))
		throwError(sol->instance, sol, "cvtSolutionToSuccessors: Converted solution is wrong");	
}

int setSEC(double *coeffs, int *indexes, CplexData *cpx, CPXCALLBACKCONTEXTptr context, SubtoursData *subData, int iterNum, Instance *inst, int nCols, bool isBenders)
{
	int retVal = 0;
	size_t n = inst->nNodes;

	// set all coeffs to 1 at the beggining so we don't have to think about them again
	for (size_t i = 0; i < nCols; i++)
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
				indexes[nnz] = (int)xpos(next, i, n);
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

double computeSuccessorsSolCost(int *successors, Instance *inst)
{
	if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (!checkSuccessorSolution(inst, successors) != 0))
		return -1.0;

	int n = (int)inst->nNodes;
	enum EdgeWeightType ewt = inst->params.edgeWeightType ;
	bool roundFlag = inst->params.roundWeights;

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

bool checkSuccessorSolution(Instance *inst, int *successors)
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
			return false;
		indexChecked[i]++;
		i = successors[i];
		count++;
	} while (count <= n);
	
	return true;
}

int WarmStart(CplexData *cpx, int *successors)
{
	Instance *inst = cpx->inst;
	size_t n = inst->nNodes;

	if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (!checkSuccessorSolution(inst, successors) != 0))
	{
		LOG(LOG_LVL_ERROR, "WarStartSuccessors: successor solution is incorrect");
		return 1;
	}

	double *ones = malloc(n * sizeof(double));
	int *indexes = malloc(n * sizeof(int));
	for (size_t i = 0; i < n; i++)
		ones[i] = 1.0;
	for (size_t i = 0; i < n; i++)
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
