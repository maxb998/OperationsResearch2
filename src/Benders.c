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

	int *bestSuccessorsSol = malloc((n + AVX_VEC_SIZE) * sizeof(int) * 2);
	__uint128_t bestCost = sol->cost;
	if (inst->params.cplexWarmStart)
		cvtSolutionToSuccessors(sol, bestSuccessorsSol);

	int iterNum = 0;
	int errCode = 0;

	while (currentTime - startTime < tlim)
	{
		if (inst->params.cplexWarmStart || ((iterNum > 0) && inst->params.cplexPatching)) // always use warm start if the "warm-starting" solution derived from a previous iteration of benders
			if ((errCode = WarmStart(&cpx, bestSuccessorsSol) != 0))
				throwError("Benders: error on WarmStart with code %d", errCode);

		// set time limit as remainig time from starting time
		clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    	currentTime = cvtTimespec2Double(timeStruct);
		CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, tlim + startTime - currentTime);

		LOG (LOG_LVL_TRACE, "Benders[%d]: Running CPXmipopt", iterNum);
		errCode = CPXmipopt(cpx.env, cpx.lp);
		if (errCode != 0)
			throwError("Benders: CPXmipopt failed with code %d", errCode);

		LOG (LOG_LVL_TRACE, "Benders[%d]: Running CPXmipopt", iterNum);
		errCode = CPXgetx(cpx.env, cpx.lp, xstar, 0, ncols - 1);
		if (errCode != 0)
			throwError("Benders: CPXgetx failed with code %d", errCode);

		double objVal;
		errCode = CPXgetobjval(cpx.env, cpx.lp, &objVal);
		if (errCode != 0)
			throwError("Benders: CPXgetobjval failed with code %d", errCode);
		
		LOG (LOG_LVL_TRACE, "Benders[%d]: Converting CPLEX solution to successors", iterNum);
		cvtCPXtoSuccessors(xstar, ncols, inst, &sub);

		if (sub.subtoursCount == 1) // means that there is only one subtour
		{
			LOG(LOG_LVL_NOTICE, "Optimal Solution found at iteration %d with cost %lf", iterNum, cvtCost2Double(computeSuccessorsSolCost(sub.successors, inst)));
			swapElems(bestSuccessorsSol, sub.successors)
			sub.subtoursMap = &sub.successors[n + AVX_VEC_SIZE];
			break;
		}

		// add subtour elimination constraints
		LOG (LOG_LVL_TRACE, "Benders[%d]: Setting subtour elimination constraints", iterNum);
		errCode = setSEC(xstar, indexes, &cpx, NULL, &sub, iterNum, inst, ncols);
		if (errCode != 0)
			throwError("Benders: SetSEC failed with code %d", errCode);
		
		if (inst->params.cplexPatching)
		{
			// generate a solution using Repair Heuristic and check if it is better than the previous solutions
			__uint128_t cost = PatchingHeuristic(&sub, inst);

			LOG(LOG_LVL_DEBUG, "Subtours at iteration %d is %d. ObjValue = %lf Cost of Repaired Solution: %lf", iterNum, sub.subtoursCount, objVal, cvtCost2Double(cost));

			if (cost < bestCost)
			{
				swapElems(bestSuccessorsSol, sub.successors)
				sub.subtoursMap = &sub.successors[n + AVX_VEC_SIZE];
				bestCost = cost;
				LOG(LOG_LVL_INFO, "Found new best solution with cost %lf", cvtCost2Double(bestCost));
			}
		}
		else
			LOG(LOG_LVL_INFO, "Subtours at iteration %d is %d. ObjValue = %lf", iterNum, sub.subtoursCount, objVal);

		iterNum++;
	}

	cvtSuccessorsToSolution(bestSuccessorsSol, bestCost, sol);

	free(xstar);
	destroyCplexData(&cpx);
	destroySubtoursData(&sub);
	free(bestSuccessorsSol);

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    currentTime = cvtTimespec2Double(timeStruct);
	sol->execTime += currentTime - startTime;
}

