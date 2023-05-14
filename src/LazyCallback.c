#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


// Callback for adding the lazy subtour elimination cuts
static int CPXPUBLIC subtourEliminationCallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) ;

Solution lazyCallback(Instance *inst)
{
	// Stores the solution so that we can populate cbData and return a Solution type for the method
	Solution sol = newSolution(inst);
	CplexData cpx = initCplexData(inst);

	int ncols = CPXgetnumcols(cpx.env, cpx.lp);

	CallbackData cbData = 
	{
		.sol = &sol,
		.cpx = &cpx,
		.ncols = ncols
	};
	
	if (CPXsetintparam(cpx.env, CPX_PARAM_MIPCBREDLP, CPX_OFF))
		cplexError(&cpx, inst, &sol, "lazyCallback: error on CPXsetinitparam(CPX_PARAM_MIPCBREDLP)");
	if (CPXcallbacksetfunc(cpx.env, cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, subtourEliminationCallback, &cbData))
		cplexError(&cpx, inst, &sol, "lazyCallback: error on CPXsetlazyconstraintcallbackfunc");
	
	/*int ncores = 1;
	if (CPXgetnumcores(cpx.env, &ncores))
		cplexError(&cpx, inst, &tempSolution, "lazyCallback: error on CPXgetnumcores");
	if (ncores != inst->params.nThreads)*/
	if (CPXsetintparam(cpx.env, CPX_PARAM_THREADS, inst->params.nThreads))
		cplexError(&cpx, inst, &sol, "lazyCallback: error on CPXsetintparam(CPX_PARAM_THREADS)");

	if (CPXmipopt(cpx.env, cpx.lp))
		cplexError(&cpx, inst, NULL, "lazyCallback: output of CPXmipopt != 0");
	
	return sol;
}

static int CPXPUBLIC subtourEliminationCallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) 
{
	CallbackData *cbData = (CallbackData *) userhandle;

	Solution *sol = cbData->sol;
	Instance *inst = sol->instance;

	// Stores the number of nodes in the problem
	int nNodes = inst->nNodes;

	// Stores current vector x from CPLEX
	double *xstar = malloc(cbData->ncols * sizeof(double));
	if (xstar == NULL) return callbackError("subtourEliminationCallback: xstar allocation error: out of memory");
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, cbData->ncols-1, NULL) != 0) return callbackError("subtourEliminationCallback: CPXgetcallbacknodex error.");

	// Stores the solution with the "successor" convention
	SubtoursData subData = {0};
	subData.successors = malloc(nNodes * 2 * sizeof(int)); 
	if (subData.successors == NULL) return callbackError("subtourEliminationCallback: successor allocation error.");
	
	// Stores the subtour in which the nodes represented buy the indexes are in
	subData.subtoursMap = &subData.successors[nNodes];

	// Stores the number of subtours foúnd by CPLEX
	subData.subtoursCount = 0;
	
	// Stores the index of the coefficients that we pass to CPLEX 
	int *indexes = malloc(cbData->ncols * sizeof(int));
	if (indexes == NULL) return callbackError("subtourEliminationCallback: indexes allocation error.");

	cvtCPXtoSuccessors(xstar, cbData->ncols, nNodes, &subData);

	if(subData.subtoursCount != 1)
	{
		double * coeffs = xstar;
		setSEC(coeffs, indexes, NULL, context, &subData, 0, inst, cbData->ncols, 0);
		if (CPXcallbackrejectcandidate(context, 0, 0, NULL, NULL, NULL, NULL, NULL))
			return callbackError("subtourEliminationCallback: candidate solution rejection failed");
	}
	else
		cvtSuccessorsToSolution(subData.successors, cbData->sol);
	

	free(xstar);
	free(subData.successors);
	free(indexes);
	return 0;
}
