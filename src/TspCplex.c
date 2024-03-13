#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <stdio.h>
#include <stdarg.h> // used for logger va_list


CplexData initCplexData(Instance *inst)
{
	CplexData cpxData = { .inst = inst };

	int errCode = 0;
	cpxData.env = CPXopenCPLEX(&errCode);
	if (errCode)
		throwError("buildCPXModel: error at CPXopenCPLEX with code %d", errCode);

	// screen output
	// if (inst->params.logLevel >= LOG_LVL_TRACE)
	// 	CPXsetintparam(cpxData.env, CPX_PARAM_SCRIND, CPX_ON);

	// random seed for cplex
	if (inst->params.randomSeed != -1)
		CPXsetintparam(cpxData.env, CPX_PARAM_RANDOMSEED, inst->params.randomSeed);
	else
		CPXsetintparam(cpxData.env, CPX_PARAM_RANDOMSEED, time(NULL));

	cpxData.lp = CPXcreateprob(cpxData.env, &errCode, "TSP");
	if (errCode)
		throwError("buildCPXModel: error at CPXcreateprob with code %d", errCode);

    int n = inst->nNodes;

	char binary = CPX_BINARY; 
	char cname[30];
	char *cnamePtr = (char*)&cname;

    // add binary var.s x(i,j) for i < j  
	
	for ( int i = 0; i < n; i++ )
	{
		for ( int j = i+1; j < n; j++ )
		{
			sprintf((char*)cname, "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			#if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
				double obj = computeEdgeCost(inst->X[i], inst->Y[i], inst->X[j], inst->Y[j], inst);
			#elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
				double obj = inst->edgeCostMat[i * n + j];
			#endif
			double ub = 1.0;

			errCode = CPXnewcols(cpxData.env, cpxData.lp, 1, &obj, NULL, &ub, &binary, &cnamePtr);
			if (errCode != 0)
				throwError("initCplexData: CPXnewcols failed with code %d", errCode);
			#ifdef DEBUG
				int cpxColNumCurrent = CPXgetnumcols(cpxData.env, cpxData.lp);
				if ( cpxColNumCurrent-1 != xpos(i,j,n) )
					throwError("initCplexData: number of column of latest added variable(%d) does not match with xpos value(%d)", cpxColNumCurrent, xpos(i,j,n));
			#endif
		}
	} 

    // add the degree constraints 

	double *value = malloc(n * (sizeof(double) + sizeof(int)));
	int *index = (int*)&value[n];
	
	for (int i = 0; i < n; i++)
		value[i] = 1.0;
	
	for ( int i = 0; i < n; i++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf((char*)cname, "degree(%d)", i+1);
		int nnz = 0;
		for ( int j = 0; j < n; j++ )
		{
			if ( j == i ) continue;
			index[nnz] = xpos(i,j,n);
			nnz++;
		}
		int izero = 0;
		errCode = CPXaddrows(cpxData.env, cpxData.lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cnamePtr);
		if (errCode != 0)
			throwError("initCplexData: CPXaddrows() failed with code %d", errCode);
	} 

	free(value);

	return cpxData;
}

void destroyCplexData(CplexData * cpxData)
{
	CPXfreeprob(cpxData->env, &cpxData->lp);
    CPXcloseCPLEX(&cpxData->env);
}

int xpos(int i, int j, int n)
{
	#ifdef DEBUG
		if (i == j)
			throwError("xpos: i == j");
	#endif
    if (i > j)
        swapElems(i, j)

    int pos = i * n + j - ((i + 1) * (i + 2)) / 2;
    return pos;
}

void cvtCPXtoSuccessors(double *xstar, int ncols, Instance *inst, SubtoursData *sub)
{
	int n = inst->nNodes;

	// reset arrays and counter
	for (int i = 0; i < (n + AVX_VEC_SIZE) * 2; i++) // exploit memory allocation
		sub->successors[i] = -1;
	sub->subtoursCount = 0;

	for (int i = 0; i < n; i++)
	{
		int succ = i;
		for (int j = 0; j < n; j++)
		{
			if ((succ != j) && (xstar[xpos(succ, j, n)] > 0.5) && (sub->subtoursMap[j] == -1))
			{
				sub->successors[succ] = j;
				sub->subtoursMap[succ] = sub->subtoursCount;
				succ = j;
				j = 0;
			}
		}
		if (succ != i)
		{
			sub->successors[succ] = i;
			sub->subtoursMap[succ] = sub->subtoursCount;
			(sub->subtoursCount)++;
		}
	}

	#ifdef DEBUG
		if (!checkSubtoursData(inst, sub))
			throwError("SubtourData converted from CPLEX is incorrect");
	#endif
}

void cvtSuccessorsToSolution(int *successors, __uint128_t cost, Solution *sol)
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

	sol->cost = cost;

	#ifdef DEBUG
		if (!checkSolution(sol))
			throwError("cvtSuccessorsToSolution: Converted solution is not correct");
	#endif
}

void cvtSolutionToSuccessors(Solution *sol, int* successors)
{
	int n = sol->instance->nNodes;

	for (int i = 0; i < n; i++)
		successors[sol->indexPath[i]] = sol->indexPath[i+1];

	#ifdef DEBUG
		if (!checkSuccessorSolution(sol->instance, successors) != 0)
			throwError("cvtSolutionToSuccessors: Converted solution is wrong");	
	#endif
}

int setSEC(double *coeffs, int *indexes, CplexData *cpx, CPXCALLBACKCONTEXTptr context, SubtoursData *subData, int iterNum, Instance *inst, int nCols)
{
	int n = inst->nNodes;

	// set all coeffs to 1 at the beggining so we don't have to think about them again
	for (int i = 0; i < nCols; i++)
		coeffs[i] = 1.;

	char sense[2] = { 'L', 0};
	char cname[30]; // size of 30 should be more than enough
	char *cnamePtr = cname;
	int izero = 0;

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
				indexes[nnz] = xpos(next, i, n);
				nnz++;
			}
			rhs += 1;
			next = subData->successors[next];
		} while (next != subtourStart);

		if (inst->params.mode & MODE_BENDERS)
			sprintf(cname, "SEC(%d,%d)", iterNum, subtourID);

		int errCode;
		if(cpx)
			errCode = CPXaddrows(cpx->env, cpx->lp, 0, 1, nnz, &rhs, sense, &izero, indexes, coeffs, NULL, &cnamePtr);
		else
			errCode = CPXcallbackrejectcandidate(context, 1, nnz, &rhs, sense, &izero, indexes, coeffs);
		
		if (errCode != 0)
			return errCode;
	}

	return 0;
}

__uint128_t computeSuccessorsSolCost(int *successors, Instance *inst)
{
	if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (!checkSuccessorSolution(inst, successors) != 0))
		throwError("computeSuccessorsSolCost: successors array passed as input is not feasible");

	int n = inst->nNodes;

	__uint128_t cost = 0LL;

	int i = 0;
	int counter = 0;
	do
	{
		int succ = successors[i];
		#if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
			cost += cvtFloat2Cost(computeEdgeCost(inst->X[i], inst->Y[i], inst->X[succ], inst->Y[succ], inst));
		#elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
			cost += cvtFloat2Cost(inst->edgeCostMat[i * (size_t)n + succ]);
		#endif
		i = succ;

		if (counter > n)
			throwError("computeSuccessorsSolCost: There are subtours inside the successor array even after repair heuristic");
		counter++;

	} while (i != 0);
	
	return cost;
}

bool checkSubtoursData(Instance *inst, SubtoursData *sub)
{
	int n = inst->nNodes;

	bool *appearedOnce = malloc(n * sizeof(bool));
	for (int i = 0; i < n; i++)
		appearedOnce[i] = false;
	for (int i = 0; i < n; i++)
	{
		int succ = sub->successors[i];
		if (succ < 0 || succ > n)
		{
			LOG(LOG_LVL_ERROR, "checkSubtoursData: successors[%d]=%d is invalid", i, succ);
			return false;
		}
		if (appearedOnce[succ])
		{
			LOG(LOG_LVL_ERROR, "checkSubtoursData: successors[%d]=%d is present more than once in the successors array", i, succ);
			return false;
		}
		appearedOnce[succ] = true;
	}
	free(appearedOnce);

	for (int i = 0; i < sub->subtoursCount; i++)
	{
		// find first element of subtour i
		int firstElem;
		for (firstElem = 0; firstElem < n; firstElem++)
			if (sub->subtoursMap[firstElem] == i)
				break;

		// check if all members of subtour i are actually mapped as part of subtour i
		for (int j = sub->successors[firstElem]; j != firstElem; j = sub->successors[j])
		{
			if (sub->subtoursMap[j] != i)
			{
				LOG(LOG_LVL_ERROR, "checkSubtoursData: element at position %d is not part of the subtour %d as it should(subtoursMap incoherent)", j, i);
				return false;
			}
		}
	}
	
	return true;
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
		{ LOG(LOG_LVL_ERROR, "checkSuccessorSolution: Node %d appears multiple times in successors", i); retVal = false; break; }

		indexChecked[i]++;
		i = successors[i];
		count++;
	} while (count <= n);

	free(indexChecked);
	
	return retVal;
}

int WarmStart(CplexData *cpx, int *successors)
{
	Instance *inst = cpx->inst;
	int n = inst->nNodes;

	#ifdef DEBUG
		if (!checkSuccessorSolution(inst, successors))
			throwError("WarStartSuccessors: successor solution is incorrect");
	#endif

	double *ones = malloc(n * (sizeof(double) + sizeof(int)));
	int *indexes = (int*)&ones[n];
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

	return retVal;
}


SubtoursData initSubtoursData(int n)
{
	SubtoursData subData = {
		.subtoursCount = 0,
		.successors = malloc((n + AVX_VEC_SIZE) * sizeof(int) * 2),
	};

	if (subData.successors == NULL)
		throwError("initSubtoursData: Failed to allocate memory");
	
	subData.subtoursMap = &subData.successors[n + AVX_VEC_SIZE];

	return subData;
}


void destroySubtoursData(SubtoursData *subData)
{
	free(subData->successors);
	// free(subData->subtoursMap);
	subData->successors = NULL;
	// subData->subtoursMap = NULL;
}
