#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// Callback for adding the lazy subtour elimination cuts
static int CPXPUBLIC subtourEliminationCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

Solution lazyCallback(Instance *inst, double tLimSec)
{
	Solution sol = newSolution(inst);
	CplexData cpx = initCplexData(inst);

	sol.ncols = CPXgetnumcols(cpx.env, cpx.lp);

	CPXsetintparam(cpx.env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	
	CPXsetlazyconstraintcallbackfunc(cpx.env, subtourEliminationCallback, &sol);
	int ncores = 1;
	CPXgetnumcores(cpx.env, &ncores);
	CPXsetintparam(cpx.env, CPX_PARAM_THREADS, ncores); 


	destroyCplexData(&cpx);
	return sol;
}

static int CPXPUBLIC subtourEliminationCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	Solution *sol = (Solution *) cbhandle;

	// Stores current vector x from CPLEX
	double *xstar = (double*) malloc(sol->ncols * sizeof(double));
	if(xstar == NULL) throwError(sol->instance, sol, "subtourEliminationCallback: xstar allocation error.");
	if(CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, sol->ncols-1) != 0) throwError(sol->instance, sol, "subtourEliminationCallback: CPXgetcallbacknodex error.");

	// Stores the solution with the "successor" convention
	int *successors = calloc(sol->instance->nNodes, sizeof(int)); 
	if(successors == NULL) throwError(sol->instance, sol, "subtourEliminationCallback: successor allocation error.");
	
	// Stores the subtour in which the nodes represented buy the indexes are in
	int *subtoursMap = calloc(sol->instance->nNodes, sizeof(int));
	if(subtoursMap == NULL) throwError(sol->instance, sol, "subtourEliminationCallback: subtoursMap allocation error.");

	// Stores ???
	int *indexes = calloc(sol->ncols, sizeof(int));
	if(indexes == NULL) throwError(sol->instance, sol, "subtourEliminationCallback: indexes alloation error.");

	// Stores the number of subtours foúnd by CPLEX
	int subtourCount = 0;
	for (size_t i = 0; i < sol->instance->nNodes; i++)
	{
		size_t succ = i;
		for (size_t j = 0; j < sol->instance->nNodes; j++)
		{
			if ((succ != j) && (xstar[xpos(succ, j, sol->instance->nNodes)] > 0.5))
			{
				successors[succ] = (int)j;
				subtoursMap[succ] = subtourCount;
				LOG(LOG_LVL_EVERYTHING, "x(%3d,%3d) = 1   subtour n° %d\n", succ, j, subtoursMap[succ]);
				succ = j;
				j = 0;
			}
		}
		if (succ != i)
		{
			successors[succ] = i;
			subtoursMap[succ] = subtourCount;
			LOG(LOG_LVL_EVERYTHING, "x(%3d,%3d) = 1   subtour n° %d\n", succ, i, subtoursMap[succ]);
			subtourCount++;
		}
	}

	free(xstar);
	free(successors);
	free(subtoursMap);
	free(indexes);
}
