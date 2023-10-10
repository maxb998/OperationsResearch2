#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time



void benders(Solution *sol, double tlim)
{
	struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);
	double currentTime = startTime;

	Instance *inst = sol->instance;
	int n = inst->nNodes;

	if (!checkSolution(sol))
		throwError("benders: Input solution is not valid");

    CplexData cpx = initCplexData(inst);

	int ncols = CPXgetnumcols(cpx.env, cpx.lp);
	double *xstar = malloc(ncols * sizeof(double));

	SubtoursData sub = initSubtoursData(n);

	int *bestSuccessorsSol = malloc(n * sizeof(int));
	__uint128_t bestCost = sol->cost;
	cvtSolutionToSuccessors(sol, bestSuccessorsSol);

	int *indexes = malloc(ncols * sizeof(int));

	int iterNum = 0;

	while (currentTime - startTime < tlim)
	{
		WarmStart(&cpx, bestSuccessorsSol);

		// set time limit as remainig time from starting time
		clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    	currentTime = cvtTimespec2Double(timeStruct);
		CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, currentTime + tlim - startTime);


		if (CPXmipopt(cpx.env, cpx.lp))
			throwError("Benders: output of CPXmipopt != 0");

		if (CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1))
			throwError("Benders: output of CPXgetx != 0");

		sub.subtoursCount = 0;
		
		cvtCPXtoSuccessors(xstar, ncols, n, &sub);

		if (sub.subtoursCount == 1) // means that there is only one subtour
		{
			LOG(LOG_LVL_LOG, "Optimal Solution found.");
			register int *temp;
			swapElems(bestSuccessorsSol, sub.successors, temp);
			break;
		}

		// add subtour elimination constraints
		double * coeffs = xstar; // reuse xStar instead of allocating new memory

		if (setSEC(coeffs, indexes, &cpx, NULL, &sub, iterNum, inst, ncols, 1) == 1)
			throwError("Benders: SetSEC failed");
		
		// generate a solution using Repair Heuristic and check if it is better than the previous solutions
		__uint128_t cost = PatchingHeuristic(&sub, inst);

		if (!checkSuccessorSolution(inst, sub.successors))
			throwError("Benders: Successors after repair heuristic does not represent a loop");

		if (cost < bestCost)
		{
			register int *temp;
			swapElems(bestSuccessorsSol, sub.successors, temp);
			bestCost = cost;
		}

		LOG(LOG_LVL_LOG, "Subtours at iteration %d is %d. Cost of Repaired Solution: %lf  Incumbent: %lf", iterNum, sub.subtoursCount, cvtCost2Double(cost), cvtCost2Double(bestCost));

		iterNum++;
	}

	free(xstar);
	free(indexes);
	destroyCplexData(&cpx);
	destroySubtoursData(&sub);

	cvtSuccessorsToSolution(bestSuccessorsSol, sol);

	free(bestSuccessorsSol);

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    currentTime = cvtTimespec2Double(timeStruct);
	sol->execTime += currentTime - startTime;
}

