#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time



void benders(Solution *sol, double tlim)
{
	struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	double currentTimeSec = startTimeSec;

	Instance *inst = sol->instance;
	size_t n = inst->nNodes;

	if (!checkSolution(sol))
		throwError(inst, sol, "benders: Input solution is not valid");

    CplexData cpx = initCplexData(inst);

	int ncols = CPXgetnumcols(cpx.env, cpx.lp);
	double *xstar = malloc(ncols * sizeof(double));

	SubtoursData sub = {
        .subtoursCount = 0,
        .successors = malloc(n * sizeof(int)),
        .subtoursMap = malloc(n * sizeof(int))
    };

	int *bestSuccessorsSol = malloc(n * sizeof(int));
	double bestCost = sol->cost;
	cvtSolutionToSuccessors(sol, bestSuccessorsSol);

	int *indexes = malloc(ncols * sizeof(int));

	int iterNum = 0;

	while (currentTimeSec - startTimeSec < tlim)
	{
		WarmStart(&cpx, bestSuccessorsSol);

		// set time limit as remainig time from starting time
		clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    	currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
		CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, currentTimeSec + tlim - startTimeSec);


		if (CPXmipopt(cpx.env, cpx.lp))
		{
			destroyCplexData(&cpx);
			throwError(inst, NULL, "Benders: output of CPXmipopt != 0");
		}

		if (CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1))
		{
			destroyCplexData(&cpx);
			throwError(inst, NULL, "Benders: output of CPXgetx != 0");
		}

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
		{
			free(xstar); free(sub.subtoursMap); free(sub.successors); free(bestSuccessorsSol); free(indexes);
			destroyCplexData(&cpx);
			throwError(inst, NULL, "Benders: SetSEC failed");
		}
		
		// generate a solution using Repair Heuristic and check if it is better than the previous solutions
		double cost = PatchingHeuristic(&sub, inst);

		if (!checkSuccessorSolution(inst, sub.successors))
		{
			free(xstar); free(sub.subtoursMap); free(sub.successors); free(bestSuccessorsSol); free(indexes);
			destroyCplexData(&cpx);
			throwError(inst, NULL, "Benders: Successors after repair heuristic does not represent a loop");
		}

		if (cost < bestCost)
		{
			register int *temp;
			swapElems(bestSuccessorsSol, sub.successors, temp);
			bestCost = cost;
		}

		LOG(LOG_LVL_LOG, "Subtours at iteration %d is %d. Cost of Repaired Solution: %lf  Incumbent: %lf", iterNum, sub.subtoursCount, cost, bestCost);

		iterNum++;
	}

	free(xstar);
	free(indexes);
	destroyCplexData(&cpx);

	free(sub.successors);
	free(sub.subtoursMap);

	cvtSuccessorsToSolution(bestSuccessorsSol, sol);

	free(bestSuccessorsSol);

	sol->cost = bestCost;
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	sol->execTime += currentTimeSec - startTimeSec;
}

