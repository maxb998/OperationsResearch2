#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


// Function to post the solution to cplex. Vals and indexes are two allocated portions of memory of ncols elements each.
static int PostSolution(CPXCALLBACKCONTEXTptr context, Instance *inst, int ncols, int *successors, __uint128_t cost, double *vals, int *indexes);

// Callback for adding the lazy subtour elimination cuts



void BranchAndCut(Solution *sol, double tlim)
{
	struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

	Instance *inst = sol->instance;

	if (!checkSolution(sol))
		throwError("benders: Input solution is not valid");

	CplexData cpx = initCplexData(inst);
	int errCode = 0;

	CallbackData cbData = initCallbackData(&cpx, sol);

	if ((errCode = WarmStart(&cpx, cbData.bestSuccessors)) != 0)
		throwError("BranchAndCut: error on WarmStart with code %d", errCode);

	if (CPXsetintparam(cpx.env, CPX_PARAM_MIPCBREDLP, CPX_OFF))
		throwError("BranchAndCut: error on CPXsetinitparam(CPX_PARAM_MIPCBREDLP)");

	if (CPXcallbacksetfunc(cpx.env, cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericCallbackCandidate, &cbData))
		throwError("BranchAndCut: error on CPXsetlazyconstraintcallbackfunc");
	
	if (CPXsetintparam(cpx.env, CPX_PARAM_THREADS, inst->params.nThreads))
		throwError("BranchAndCut: error on CPXsetintparam(CPX_PARAM_THREADS)");

	if (CPXmipopt(cpx.env, cpx.lp))
		throwError("BranchAndCut: output of CPXmipopt != 0");
	
	cvtSuccessorsToSolution(cbData.bestSuccessors, sol);
	//sol->cost = cbData.bestCost;

	destroyCallbackData(&cbData);

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
	sol->execTime += cvtTimespec2Double(timeStruct) - startTime;
}

int CPXPUBLIC genericCallbackCandidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) 
{
	CallbackData *cbData = (CallbackData *) userhandle;

	Instance *inst = cbData->inst;
	int errCode = 0;

	// Stores the number of nodes in the problem
	int nNodes = inst->nNodes;

	// Stores current vector x from CPLEX
	double *xstar = malloc(cbData->ncols * sizeof(double));
	if ((errCode = CPXcallbackgetcandidatepoint(context, xstar, 0, cbData->ncols-1, NULL)) != 0) 
	{
		LOG(LOG_LVL_ERROR, "subtourEliminationCallback: CPXgetcallbacknodex failed with code %d", errCode);
		return 1;
	}

	SubtoursData sub = initSubtoursData(inst->nNodes);
	
	// Stores the index of the coefficients that we pass to CPLEX 
	int *indexes = malloc(cbData->ncols * sizeof(int));

	cvtCPXtoSuccessors(xstar, cbData->ncols, nNodes, &sub);

	pthread_mutex_lock(&cbData->mutex);
	if (inst->params.mode == MODE_BRANCH_CUT) LOG(LOG_LVL_LOG, "Iteration %d subtours detected %d", cbData->iterNum, sub.subtoursCount);
	cbData->iterNum++;
	pthread_mutex_unlock(&cbData->mutex);

	__uint128_t cost = -1LL;
	if(sub.subtoursCount > 1)
	{
		double * coeffs = xstar;
		if ((errCode = setSEC(coeffs, indexes, NULL, context, &sub, 0, inst, cbData->ncols, 0)) != 0)
			throwError("BranchAndCut callback: setSEC failed with code %d", errCode);

		cost = PatchingHeuristic(&sub, inst);
	}
	else
		cost = computeSuccessorsSolCost(sub.successors, inst);

	// set new incumbent if found
	if (cost < cbData->bestCost)
	{
		int *temp;
		swapElems(cbData->bestSuccessors, sub.successors, temp);

		enum LogLevel lvl = LOG_LVL_NOTICE;
		if (inst->params.mode != MODE_BRANCH_CUT) lvl = LOG_LVL_LOG;
		LOG(lvl, "Branch & Cut: Incumbent updated -> cost = %lf", cvtCost2Double(cost));

		cbData->bestCost = cost;

		// post solution to cplex
		if ((sub.subtoursCount > 1) && ((errCode = PostSolution(context, inst, cbData->ncols, cbData->bestSuccessors, cbData->bestCost, xstar, indexes)) != 0))
			throwError("BranchAndCut callback: postSolution failed with code %d", errCode);
	}

	free(xstar); destroySubtoursData(&sub); free(indexes);
	return 0;
}

CallbackData initCallbackData(CplexData *cpx, Solution *sol)
{
	CallbackData cbData = {
		.inst = cpx->inst,
		.ncols = CPXgetnumcols(cpx->env, cpx->lp),
		.iterNum = 0,
		.bestCost = -1LL,
	};

	if (cbData.ncols <= 0)
		throwError("initCallbackData: CPXgetnumcols returned 0. Cpx.env is probably empty");

	cbData.bestSuccessors = malloc((cpx->inst->nNodes + AVX_VEC_SIZE) * sizeof(int));
	if (cbData.bestSuccessors == NULL)
		throwError("initCallbackData: Failed to allocate memory");

	pthread_mutex_init(&cbData.mutex, NULL);

	cvtSolutionToSuccessors(sol, cbData.bestSuccessors);
	cbData.bestCost = sol->cost;

	return cbData;
}

void destroyCallbackData(CallbackData *cbData)
{
	free(cbData->bestSuccessors);
	cbData->bestSuccessors = NULL;
	pthread_mutex_destroy(&cbData->mutex);
}

static int PostSolution(CPXCALLBACKCONTEXTptr context, Instance *inst, int ncols, int *successors, __uint128_t cost, double *vals, int *indexes)
{
	if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (!checkSuccessorSolution(inst, successors) != 0))
		return 1;

	for (int i = 0; i < ncols; i++)
		vals[i] = 0.0;
	for (int i = 0; i < ncols; i++)
		indexes[i] = i;
	for (int i = 0; i < inst->nNodes; i++)
	{
		int pos = xpos(i, successors[i], inst->nNodes);
		vals[pos] = 1.0;
	}

	int strat = CPXCALLBACKSOLUTION_NOCHECK;
	if (inst->params.logLevel >= LOG_LVL_DEBUG)
		strat = CPXCALLBACKSOLUTION_CHECKFEAS;

	int retVal = CPXcallbackpostheursoln(context, ncols, indexes, vals, cvtCost2Double(cost), strat);

	return retVal;
}
