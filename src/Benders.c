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

    CplexData cpx = initCplexData(inst);

	int ncols = CPXgetnumcols(cpx.env, cpx.lp);
	double *xstar = malloc(ncols * (sizeof(double) + sizeof(int)));
	int *indexes = (int*)&xstar[ncols];

	SubtoursData sub = initSubtoursData(n);

	#ifdef DEBUG
		if (inst->params.cplexWarmStart && (!checkSolution(sol)))
			throwError("benders: Input solution is not valid");
	#endif

	int *bestSuccessorsSol = malloc(n * sizeof(int));
	__uint128_t bestCost = sol->cost;
	if (inst->params.cplexWarmStart)
		cvtSolutionToSuccessors(sol, bestSuccessorsSol);

	int iterNum = 0;
	int errCode = 0;

	while (currentTime - startTime < tlim)
	{
		if (inst->params.cplexWarmStart || (iterNum > 0)) // always use warm start if the "warm-starting" solution derived from a previous iteration of benders
			if ((errCode = WarmStart(&cpx, bestSuccessorsSol) != 0))
				throwError("Benders: error on WarmStart with code %d", errCode);

		// set time limit as remainig time from starting time
		clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    	currentTime = cvtTimespec2Double(timeStruct);
		CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, tlim + startTime - currentTime);

		errCode = CPXmipopt(cpx.env, cpx.lp);
		if (errCode != 0)
			throwError("Benders: output of CPXmipopt != 0");

		errCode = CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1);
		if (errCode != 0)
			throwError("Benders: output of CPXgetx != 0");

		sub.subtoursCount = 0;
		
		cvtCPXtoSuccessors(xstar, ncols, n, &sub);

		if (sub.subtoursCount == 1) // means that there is only one subtour
		{
			LOG(LOG_LVL_NOTICE, "Optimal Solution found at iteration %d with cost %lf", iterNum, cvtCost2Double(computeSuccessorsSolCost(sub.successors, inst)));
			swapElems(bestSuccessorsSol, sub.successors)
			break;
		}

		// add subtour elimination constraints
		errCode = setSEC(xstar, indexes, &cpx, NULL, &sub, iterNum, inst, ncols);
		if (errCode != 0)
			throwError("Benders: SetSEC failed");
		
		// generate a solution using Repair Heuristic and check if it is better than the previous solutions
		__uint128_t cost = PatchingHeuristic(&sub, inst);

		if (!checkSuccessorSolution(inst, sub.successors))
			throwError("Benders: Successors after repair heuristic does not represent a loop");

		LOG(LOG_LVL_DEBUG, "Subtours at iteration %d is %d. Cost of Repaired Solution: %lf", iterNum, sub.subtoursCount, cvtCost2Double(cost));

		if (cost < bestCost)
		{
			swapElems(bestSuccessorsSol, sub.successors)
			bestCost = cost;
			LOG(LOG_LVL_LOG, "Found new best solution with cost %lf", cvtCost2Double(bestCost));
		}

		iterNum++;
	}

	free(xstar);
	destroyCplexData(&cpx);
	destroySubtoursData(&sub);

	cvtSuccessorsToSolution(bestSuccessorsSol, sol);

	free(bestSuccessorsSol);

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    currentTime = cvtTimespec2Double(timeStruct);
	sol->execTime += currentTime - startTime;
}

