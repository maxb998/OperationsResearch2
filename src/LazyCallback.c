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
	Solution tempSolution = newSolution(inst);
	CplexData cpx = initCplexData(inst);

	int ncols = CPXgetnumcols(cpx.env, cpx.lp);

	CallbackData cbData = 
	{
		.sol = &tempSolution,
		.cpx = &cpx,
		.ncols = ncols
	};
	
	if (CPXsetintparam(cpx.env, CPX_PARAM_MIPCBREDLP, CPX_OFF))
		cplexError(&cpx, inst, &tempSolution, "lazyCallback: error on CPXsetinitparam(CPX_PARAM_MIPCBREDLP)");
	if (CPXcallbacksetfunc(cpx.env, cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, subtourEliminationCallback, &cbData))
		cplexError(&cpx, inst, &tempSolution, "lazyCallback: error on CPXsetlazyconstraintcallbackfunc");
	
	/*int ncores = 1;
	if (CPXgetnumcores(cpx.env, &ncores))
		cplexError(&cpx, inst, &tempSolution, "lazyCallback: error on CPXgetnumcores");
	if (ncores != inst->params.nThreads)*/
	if (inst->params.logLevel >= LOG_LVL_DEBUG)
		if (CPXsetintparam(cpx.env, CPX_PARAM_THREADS, inst->params.nThreads))
			cplexError(&cpx, inst, &tempSolution, "lazyCallback: error on CPXsetintparam(CPX_PARAM_THREADS)");

	if (CPXmipopt(cpx.env, cpx.lp))
		cplexError(&cpx, inst, NULL, "lazyCallback: output of CPXmipopt != 0");

	/*double *xstar = malloc(ncols * sizeof(double));

	if (CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1))
		cplexError(&cpx, inst, NULL, "lazyCallback: output of CPXgetx != 0");*/
	


	destroyCplexData(&cpx);
	return tempSolution;
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
	if (xstar == NULL) cplexError(NULL, inst, sol, "subtourEliminationCallback: xstar allocation error: out of memory");
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, cbData->ncols-1, NULL) != 0) cplexError(NULL, inst, sol, "subtourEliminationCallback: CPXgetcallbacknodex error.");

	// Stores the solution with the "successor" convention
	int *successors = malloc(nNodes * 2 * sizeof(int)); 
	if (successors == NULL) cplexError(NULL, inst, sol, "subtourEliminationCallback: successor allocation error.");
	
	// Stores the subtour in which the nodes represented buy the indexes are in
	int *subtoursMap = &successors[nNodes];

	// Stores the number of subtours foúnd by CPLEX
	int subtourCount = 0;
	
	// Stores the index of the coefficients that we pass to CPLEX 
	int *indexes = malloc(cbData->ncols * sizeof(int));
	if (indexes == NULL) cplexError(NULL, inst, sol, "subtourEliminationCallback: indexes allocation error.");

	cvtCPXtoSuccessors(xstar, cbData->ncols, nNodes, successors, subtoursMap, &subtourCount);

	if(subtourCount != 1)
	{
		double * coeffs = xstar;
		setSEC(coeffs, indexes, cbData->cpx, successors, subtoursMap, subtourCount, 0, inst, cbData->ncols, 0);
		if (CPXcallbackrejectcandidate(context, 0, 0, NULL, NULL, NULL, NULL, NULL))
			cplexError(NULL, inst, sol, "subtourEliminationCallback: candidate solution rejection failed");
	}
	

	free(xstar);
	free(successors);
	free(indexes);
	return 0;
}
