#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// Callback for adding the lazy subtour elimination cuts
static int CPXPUBLIC subtourEliminationCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

Solution lazyCallback(Instance *inst, double tLimSec)
{
	// Stores the solution so that we can populate cbData and return a Solution type for the method
	Solution tempSolution = newSolution(inst);
	CplexData cpx = initCplexData(inst);

	CallbackData cbData = 
	{
		.sol = &tempSolution,
		.cpx = &cpx,
		.ncols = CPXgetnumcols(cpx.env, cpx.lp)
	};
	

	CPXsetintparam(cpx.env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	
	CPXsetlazyconstraintcallbackfunc(cpx.env, subtourEliminationCallback, &cbData);
	int ncores = 1;
	CPXgetnumcores(cpx.env, &ncores);
	CPXsetintparam(cpx.env, CPX_PARAM_THREADS, ncores); 


	destroyCplexData(&cpx);
	return tempSolution;
}

static int CPXPUBLIC subtourEliminationCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	CallbackData *cbData = (CallbackData *) cbhandle;

	// Stores the number of nodes in the problem
	int nNodes = cbData->sol->instance->nNodes;

	// Stores current vector x from CPLEX
	double *xstar = (double*) malloc(cbData->ncols * sizeof(double));
	if(xstar == NULL) throwError(cbData->sol->instance, &cbData->sol, "subtourEliminationCallback: xstar allocation error.");
	if(CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, cbData->ncols-1) != 0) throwError(cbData->sol->instance, &cbData->sol, "subtourEliminationCallback: CPXgetcallbacknodex error.");

	// Stores the solution with the "successor" convention
	int *successors = malloc(nNodes * sizeof(int)); 
	if(successors == NULL) throwError(cbData->sol->instance, &cbData->sol, "subtourEliminationCallback: successor allocation error.");
	
	// Stores the subtour in which the nodes represented buy the indexes are in
	int *subtoursMap = malloc(nNodes * sizeof(int));
	if(subtoursMap == NULL) throwError(cbData->sol->instance, &cbData->sol, "subtourEliminationCallback: subtoursMap allocation error.");

	// Stores the number of subtours foúnd by CPLEX
	int subtourCount = 0;
	
	// Stores the index of the coefficients that we pass to CPLEX 
	int *indexes = calloc(cbData->ncols, sizeof(int));
	if(indexes == NULL) throwError(cbData->sol->instance, &cbData->sol, "subtourEliminationCallback: indexes alloation error.");

	cvtCPXtoSuccessors(xstar, cbData->ncols, nNodes, successors, subtoursMap, subtourCount);

	if(subtourCount != 1)
	{
		double * coeffs = xstar;
		setSEC(coeffs, indexes, cbData->cpx, successors, subtoursMap, subtourCount, 0, cbData->sol->instance, cbData->ncols, 0);
	}
	

	free(xstar);
	free(successors);
	free(subtoursMap);
	free(indexes);
}
