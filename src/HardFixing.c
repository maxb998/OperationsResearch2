#include "TspCplex.h"
#include "Tsp.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// Reset Upper and Lower bounds of the problem in cplex. Pass already allocated memory in indexes and bounds (ncols elements each)
static int resetBounds(CplexData *cpx, int ncols, int *indexes, char *boundsType, double *bounds);

// Fix variables in random way up the the fixAmount, but not all the way for sure.
static int randomFix(CplexData *cpx, int fixAmount, int *successors, int *indexes, char *boundsType, double *bounds);


void HardFixing(Solution *sol, double fixingAmount, enum HardFixPolicy policy, double tlim)
{
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	double currentTimeSec = startTimeSec;

    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    int errCode = 0;

    int fixAmount = (int)((double)n * fixingAmount);

    if (inst->params.logLevel >= LOG_LVL_DEBUG)
        checkSolution(sol);

    CplexData cpx = initCplexData(inst);
    int ncols = CPXgetnumcols(cpx.env, cpx.lp);

    // allocate necessary memory
    double *xstar = malloc(ncols * sizeof(double));
    int *indexes = malloc(ncols * sizeof(int));
    char *boundsType = malloc(ncols);
    SubtoursData sub = {
        .subtoursCount = 0,
        .successors = malloc(n * sizeof(int)),
        .subtoursMap = malloc(n * sizeof(int))
    };
    CallbackData cbData = {
        .bestCost = INFINITY,
        .bestSuccessors = malloc(n * sizeof(int)),
        .ncols = ncols,
        .inst = inst,
        .iterNum = 0
    };
    cvtSolutionToSuccessors(sol, cbData.bestSuccessors);
    cbData.bestCost = sol->cost;

    if ((errCode = CPXcallbacksetfunc(cpx.env, cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericCallbackCandidate, &cbData)) != 0)
	{
		destroyCplexData(&cpx); free(xstar); free(indexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
		throwError(inst, NULL, "HardFix: error on CPXsetlazyconstraintcallbackfunc with code %d", errCode);
	}

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;

    while (currentTimeSec < startTimeSec + tlim)
    {
        switch (policy)
        {
        case HARDFIX_POLICY_RANDOM:
            if ((errCode = randomFix(&cpx, fixAmount, cbData.bestSuccessors, indexes, boundsType, xstar)) != 0)
            {
                destroyCplexData(&cpx); free(xstar); free(indexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
                throwError(inst, sol, "HardFix: randomFix failed with code %d", errCode);
            }
            break;
        case HARDFIX_POLICY_SMALLEST:
            break;
        case HARDFIX_POLICY_PROBABILITY:
            break;
        case HARDFIX_POLICY_MIXED:
            break;
        }

        double remainingTime = startTimeSec + tlim - currentTimeSec;
        if ((errCode = CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, remainingTime)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: CPXsetdblparam failed with code %d", errCode);
        }

        if ((errCode = WarmStart(&cpx, cbData.bestSuccessors)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: CPXsetdblparam failed with code %d", errCode);
        }

        if ((errCode = CPXmipopt(cpx.env, cpx.lp)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: CPXsetdblparam failed with code %d", errCode);
        }

        if ((errCode = resetBounds(&cpx, ncols, indexes, boundsType, xstar)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: resetBounds failed with code %d", errCode);
        }

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    }
    destroyCplexData(&cpx);
    free(xstar); free(indexes); free(boundsType); free(sub.successors); free(sub.subtoursMap);
    
    if (cbData.bestCost < sol->cost)
    {
        cvtSuccessorsToSolution(cbData.bestSuccessors, sol);
        sol->cost = computeSolutionCost(sol);
    }
    else
        LOG(LOG_LVL_WARNING, "HardFixing: Solution could not be optimized any further");

    free(cbData.bestSuccessors);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    sol->execTime += currentTimeSec - startTimeSec;
}

static int resetBounds(CplexData *cpx, int ncols, int *indexes, char *boundsType, double *bounds)
{
    // set lower bounds
    for (int i = 0; i < ncols; i++)
    {
        boundsType[i] = 'L';
        indexes[i] = i;
        bounds[i] = 0.0;
    }
    int retVal = CPXchgbds(cpx->env, cpx->lp, ncols, indexes, boundsType, bounds);
    if (retVal != 0)  return retVal;

    // set upper bounds
    for (int i = 0; i < ncols; i++)
    {
        boundsType[i] = 'U';
        bounds[i] = 1.0;
    }
    retVal = CPXchgbds(cpx->env, cpx->lp, ncols, indexes, boundsType, bounds);

    return retVal;
}

static int randomFix(CplexData *cpx, int fixAmount, int *successors, int *indexes, char *boundsType, double *bounds)
{
    size_t n = cpx->inst->nNodes;

    int threshold = (int)((double)fixAmount/(double)n * (double)RAND_MAX);

    // set lower bounds
    int i = 0;
    for (int j = 0; j < n && i < fixAmount; j++)
    {
        int rnd = rand();
        if (rnd < threshold)
        {
            boundsType[i] = 'B';
            indexes[i] = xpos(j, successors[j], n);
            bounds[i] = 1.0;
            i++;
        }
    }
    int retVal = CPXchgbds(cpx->env, cpx->lp, i, indexes, boundsType, bounds);
    LOG(LOG_LVL_DEBUG, "The number of fixed elements is %d and the fixAmount is set to %d", i, fixAmount);

    return retVal;
}

