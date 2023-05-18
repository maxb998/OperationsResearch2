#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time



Solution benders(Instance *inst, double tlim)
{
	size_t n = inst->nNodes;

    CplexData cpx = initCplexData(inst);

	struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	double currentTimeSec = startTimeSec;

	int ncols = CPXgetnumcols(cpx.env, cpx.lp);
	double *xstar = malloc(ncols * sizeof(double));

	SubtoursData subData = {0};
	subData.successors = malloc(n * sizeof(int));
	subData.subtoursMap = malloc(n * sizeof(int));

	double bestCost = INFINITY;
	int *bestSuccessorsSol = malloc(n * sizeof(int));

	int *indexes = malloc(ncols * sizeof(int));

	while (currentTimeSec - startTimeSec < tlim)
	{
		static int iterNum = 0;

		// set time limit as remainig time from starting time
		clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    	currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
		CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, currentTimeSec + tlim - startTimeSec);

		if (CPXmipopt(cpx.env, cpx.lp))
			cplexError(&cpx, inst, NULL, "Benders: output of CPXmipopt != 0");

		if (CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1))
			cplexError(&cpx, inst, NULL, "Benders: output of CPXgetx != 0");

		subData.subtoursCount = 0;
		
		cvtCPXtoSuccessors(xstar, ncols, n, &subData);

		if (subData.subtoursCount == 1) // means that there is only one subtour
		{
			LOG(LOG_LVL_LOG, "Optimal Solution found.");
			register int *temp;
			swapElems(bestSuccessorsSol, subData.successors, temp);
			break;
		}

		// add subtour elimination constraints
		double * coeffs = xstar; // reuse xStar instead of allocating new memory

		setSEC(coeffs, indexes, &cpx, NULL, &subData, iterNum, inst, ncols, 1);
		
		// generate a solution using Repair Heuristic and check if it is better than the previous solutions
		double cost = RepairHeuristicSuccessors(&subData, inst);

		if (checkSuccessorSolution(inst, subData.successors))
		{
			free(xstar); free(subData.subtoursMap); free(subData.successors); free(bestSuccessorsSol); free(indexes);
			cplexError(&cpx, inst, NULL, "Benders: Successors after repair heuristic does not represent a loop");
		}

		if (cost < bestCost)
		{
			register int *temp;
			swapElems(bestSuccessorsSol, subData.successors, temp);
			bestCost = cost;
		}

		LOG(LOG_LVL_LOG, "Subtours at iteration %d is %d. Cost of Repaired Solution: %lf  Incumbent: %lf", iterNum, subData.subtoursCount, cost, bestCost);

		iterNum++;
	}

	free(xstar);
	free(indexes);
	destroyCplexData(&cpx);

	free(subData.successors);
	free(subData.subtoursMap);

	Solution sol = newSolution(inst);

	cvtSuccessorsToSolution(bestSuccessorsSol, &sol);

	free(bestSuccessorsSol);

	sol.cost = computeSolutionCost(&sol);
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	sol.execTime = currentTimeSec - startTimeSec;

    return sol;
}

