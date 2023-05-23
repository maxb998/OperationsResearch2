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

	double bestCost;
	int *bestSuccessors;
} CallbackData;



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
		.bestSuccessors = malloc(n * sizeof(int)),
		.bestCost = INFINITY,
	};

	if (warmStartSol)
		WarmStart(&cpx, warmStartSol);

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
		cbData->bestCost = cost;
		int *temp;
		swapElems(cbData->bestSuccessors, sub.successors, temp);

		// post solution to cplex
		if ((sub.subtoursCount > 1) && ((errCode = PostSolution(context, inst, cbData->bestSuccessors, cbData->bestCost)) != 0))
		{
			free(xstar); free(sub.successors); free(sub.subtoursMap); free(indexes);
			LOG(LOG_LVL_ERROR,"BranchAndCut callback: postSolution failed with code %d", errCode);
			return 1;
		}

	}

	free(xstar); free(sub.successors); free(sub.subtoursMap); free(indexes);
	return 0;
}


