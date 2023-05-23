#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <pthread.h>

typedef struct
{
	Instance *inst;
	CplexData *cpx;

	int ncols; // Number of columns of the matrix inside cplex linear problem
	int iterNum;
	pthread_mutex_t mutex;

	double bestCost;
	int *bestSuccessors;
} CallbackData;

// Function to post the solution to cplex. Vals and indexes are two allocated portions of memory of ncols elements each.
static int PostSolution(CPXCALLBACKCONTEXTptr context, Instance *inst, int ncols, int *successors, double cost, double *vals, int *indexes);

// Callback for adding the lazy subtour elimination cuts
static int CPXPUBLIC genericCallbackCandidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) ;


Solution BranchAndCut(Instance *inst, double tlim, Solution *warmStartSol)
{
	struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;

	size_t n = inst->nNodes;
	CplexData cpx = initCplexData(inst);
	int ncols = CPXgetnumcols(cpx.env, cpx.lp);

	CallbackData cbData = 
	{
		.inst = inst,
		.cpx = &cpx,
		.ncols = ncols,
		.iterNum = 0,
		.bestSuccessors = malloc(n * sizeof(int)),
		.bestCost = INFINITY,
	};
	pthread_mutex_init(&cbData.mutex, NULL);

	if ((warmStartSol) && (WarmStart(&cpx, warmStartSol) != 0))
	{
		destroyCplexData(&cpx);
		throwError(inst, NULL, "lazyCallback: error on WarmStart");
	}

	if (CPXsetintparam(cpx.env, CPX_PARAM_MIPCBREDLP, CPX_OFF))
	{
		destroyCplexData(&cpx);
		throwError(inst, NULL, "lazyCallback: error on CPXsetinitparam(CPX_PARAM_MIPCBREDLP)");
	}

	if (CPXcallbacksetfunc(cpx.env, cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericCallbackCandidate, &cbData))
	{
		destroyCplexData(&cpx);
		throwError(inst, NULL, "lazyCallback: error on CPXsetlazyconstraintcallbackfunc");
	}
	
	if (CPXsetintparam(cpx.env, CPX_PARAM_THREADS, inst->params.nThreads))
	{
		destroyCplexData(&cpx);
		throwError(inst, NULL, "lazyCallback: error on CPXsetintparam(CPX_PARAM_THREADS)");
	}

	if (CPXmipopt(cpx.env, cpx.lp))
	{
		destroyCplexData(&cpx);
		throwError(inst, NULL, "lazyCallback: output of CPXmipopt != 0");
	}

	pthread_mutex_destroy(&cbData.mutex);
	
	Solution sol = newSolution(inst);
	cvtSuccessorsToSolution(cbData.bestSuccessors, &sol);
	sol.cost = cbData.bestCost;

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
	sol.execTime = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0 - startTimeSec;

	return sol;
}

static int CPXPUBLIC genericCallbackCandidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) 
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

	SubtoursData sub = {
		.successors = malloc(nNodes * sizeof(int)),
		.subtoursMap = malloc(nNodes * sizeof(int)),
		.subtoursCount = 0
	};
	
	// Stores the index of the coefficients that we pass to CPLEX 
	int *indexes = malloc(cbData->ncols * sizeof(int));

	cvtCPXtoSuccessors(xstar, cbData->ncols, nNodes, &sub);

	pthread_mutex_lock(&cbData->mutex);
	LOG(LOG_LVL_LOG, "Iteration %d subtours detected %d", cbData->iterNum, sub.subtoursCount);
	cbData->iterNum++;
	pthread_mutex_unlock(&cbData->mutex);

	double cost = INFINITY;
	if(sub.subtoursCount > 1)
	{
		double * coeffs = xstar;
		if ((errCode = setSEC(coeffs, indexes, NULL, context, &sub, 0, inst, cbData->ncols, 0)) != 0)
		{
			free(xstar); free(sub.successors); free(sub.subtoursMap); free(indexes);
			LOG(LOG_LVL_ERROR,"BranchAndCut callback: setSEC failed with code %d", errCode);
			return 1;
		}

		cost = PatchingHeuristic(&sub, inst);
	}
	else
		cost = computeSuccessorsSolCost(sub.successors, inst);

	// set new incumbent if found
	if (cost < cbData->bestCost)
	{
		int *temp;
		swapElems(cbData->bestSuccessors, sub.successors, temp);
		LOG(LOG_LVL_NOTICE, "New incumbent %lf    Old incumbent %lf", cost, cbData->bestCost);
		cbData->bestCost = cost;

		// post solution to cplex
		if ((sub.subtoursCount > 1) && ((errCode = PostSolution(context, inst, cbData->ncols, cbData->bestSuccessors, cbData->bestCost, xstar, indexes)) != 0))
		{
			free(xstar); free(sub.successors); free(sub.subtoursMap); free(indexes);
			LOG(LOG_LVL_ERROR,"BranchAndCut callback: postSolution failed with code %d", errCode);
			return 1;
		}

	}

	free(xstar); free(sub.successors); free(sub.subtoursMap); free(indexes);
	return 0;
}

static int PostSolution(CPXCALLBACKCONTEXTptr context, Instance *inst, int ncols, int *successors, double cost, double *vals, int *indexes)
{
	if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (checkSuccessorSolution(inst, successors) != 0))
		return 1;

	for (size_t i = 0; i < ncols; i++)
		vals[i] = 0.0;
	for (size_t i = 0; i < ncols; i++)
		indexes[i] = i;
	for (size_t i = 0; i < inst->nNodes; i++)
	{
		int pos = xpos(i, successors[i], inst->nNodes);
		vals[pos] = 1.0;
	}

	int strat = CPXCALLBACKSOLUTION_NOCHECK;
	if (inst->params.logLevel >= LOG_LVL_DEBUG)
		strat = CPXCALLBACKSOLUTION_CHECKFEAS;

	int retVal = CPXcallbackpostheursoln(context, ncols, indexes, vals, cost, strat);

	return retVal;
}
