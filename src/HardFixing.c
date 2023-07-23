#include "TspCplex.h"
#include "Tsp.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// Amount of nodes of the instance below which no nodes will be fixed, effectively running branch&cut. Also in any point the algorithm fixAmount won't be going any lower than this
#define MIN_UNFIX 150
// Minimum amount of edges of the solution to fix(it doesn't make much sense fixing only 1 node, if nNodes = 151 then 50 nodes will be fixed and 101 will be free)
#define MIN_FIX 50
// Maximum amount of unfixed/free edges.
#define MAX_UNFIX 300
// Incremental/Decremental step in fixAmount during computation
#define FIX_OFFSET 10
// Number of non-improving iterations before fixAmount is increased
#define STATIC_COST_THRESHOLD 10


// Reset Upper and Lower bounds of the problem in cplex. Pass already allocated memory in indexes and bounds (ncols elements each)
static int resetBounds(CplexData *cpx, int *successors, int *indexes, char *boundsType, double *bounds);

// Fix variables in random way up the the fixAmount, but not all the way for sure.
static int randomFix(CplexData *cpx, size_t fixAmount, int *fixedPositions, int *successors, int *indexes, char *boundsType, double *bounds);


void HardFixing(Solution *sol, double fixingAmount, enum HardFixPolicy policy, double tlim)
{
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	double currentTimeSec = startTimeSec;

    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    int errCode = 0;

    int fixAmount;
    if(n < MIN_UNFIX)
    {
        LOG(LOG_LVL_WARNING, "HardFix: Solution is small, so no edge will be fixed, resulting in a computation that is the same as branch-cut");
        fixAmount = 0;
    }
    else if (n < MIN_UNFIX + MIN_FIX)
        fixAmount = MIN_FIX;
    else
        fixAmount = n - MIN_UNFIX;


    if (inst->params.logLevel >= LOG_LVL_DEBUG)
        checkSolution(sol);

    CplexData cpx = initCplexData(inst);
    int ncols = CPXgetnumcols(cpx.env, cpx.lp);

    // allocate necessary memory
    double *xstar = malloc(ncols * sizeof(double));
    int *indexes = malloc(n * sizeof(int));
    int *fixesIndexes = malloc(n * sizeof(int)); for (size_t i = 0; i < n; i++) { fixesIndexes[i] = (int)i; }    
    char *boundsType = malloc(n);
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
		destroyCplexData(&cpx); free(xstar); free(indexes); free(fixesIndexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
		throwError(inst, NULL, "HardFix: error on CPXsetlazyconstraintcallbackfunc with code %d", errCode);
	}

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;

    int staticCostCount = 0;

    while (currentTimeSec < startTimeSec + tlim)
    {
        if (fixAmount > 0)
        {
            switch (policy)
            {
            case HARDFIX_POLICY_RANDOM:
                if ((errCode = randomFix(&cpx, (size_t)fixAmount, fixesIndexes, cbData.bestSuccessors, indexes, boundsType, xstar)) != 0)
                {
                    destroyCplexData(&cpx); free(xstar); free(indexes); free(fixesIndexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
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
        }

        double remainingTime = startTimeSec + tlim - currentTimeSec;
        if ((errCode = CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, remainingTime)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(fixesIndexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: CPXsetdblparam failed with code %d", errCode);
        }

        if ((errCode = WarmStart(&cpx, cbData.bestSuccessors)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(fixesIndexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: CPXsetdblparam failed with code %d", errCode);
        }

        double oldCost = cbData.bestCost;

        if ((errCode = CPXmipopt(cpx.env, cpx.lp)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(fixesIndexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: CPXsetdblparam failed with code %d", errCode);
        }

        if (oldCost <= cbData.bestCost)
            staticCostCount++;
        else
            staticCostCount = 0;

        if (fixAmount == 0)
        {
            LOG(LOG_LVL_NOTICE, "HardFix: Found the best solution");
            break;
        }
        
        if ((fixAmount > (int)n - MAX_UNFIX) && (staticCostCount > STATIC_COST_THRESHOLD))
        {
            int oldFix = fixAmount;
            fixAmount -= FIX_OFFSET;
            if (fixAmount < 0)
                fixAmount = (int)n - MAX_UNFIX;
            LOG(LOG_LVL_LOG, "Decreasing the amount of fixed edges from %lu to %lu", oldFix, fixAmount);
            staticCostCount = 0;//(n - fixAmount - MIN_UNFIX) / FIX_OFFSET;
        }

        if ((errCode = resetBounds(&cpx, cbData.bestSuccessors, indexes, boundsType, xstar)) != 0)
        {
            destroyCplexData(&cpx); free(xstar); free(indexes); free(fixesIndexes); free(boundsType); free(sub.successors); free(sub.subtoursMap); free(cbData.bestSuccessors);
            throwError(inst, sol, "HardFix: resetBounds failed with code %d", errCode);
        }

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    }
    destroyCplexData(&cpx);
    free(xstar); free(indexes); free(fixesIndexes); free(boundsType); free(sub.successors); free(sub.subtoursMap);
    
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

static int resetBounds(CplexData *cpx, int *successors, int *indexes, char *boundsType, double *bounds)
{
    size_t n = cpx->inst->nNodes;

    // set lower bounds
    for (size_t i = 0; i < n; i++)
    {
        boundsType[i] = 'L';
        indexes[i] = xpos(i, (size_t)successors[i], n);
        bounds[i] = 0.0;
    }
    int retVal = CPXchgbds(cpx->env, cpx->lp, (int)n, indexes, boundsType, bounds);
    if (retVal != 0)  return retVal;

    // set upper bounds
    for (size_t i = 0; i < n; i++)
    {
        boundsType[i] = 'U';
        bounds[i] = 1.0;
    }
    retVal = CPXchgbds(cpx->env, cpx->lp, (int)n, indexes, boundsType, bounds);

    return retVal;
}

static int randomFix(CplexData *cpx, size_t fixAmount, int *fixedPositions, int *successors, int *indexes, char *boundsType, double *bounds)
{
    size_t n = cpx->inst->nNodes;

    // n random permutations
    for (size_t i = 0; i < n; i++)
    {
        size_t pos1 = (size_t)rand() * fixAmount / RAND_MAX, pos2 = (size_t)rand() * fixAmount / RAND_MAX;
        while (pos2 == pos1) pos2 = (size_t)rand() * fixAmount / RAND_MAX;

        register int temp;
        swapElems(fixedPositions[pos1], fixedPositions[pos2], temp);
    }
    
    // setup arrays to fix bounds
    for (size_t i = 0; i < fixAmount; i++)
    {
        boundsType[i] = 'B';
        indexes[i] = xpos((size_t)fixedPositions[i], (size_t)successors[fixedPositions[i]], n);
        bounds[i] = 1.0;
    }
    int retVal = CPXchgbds(cpx->env, cpx->lp, fixAmount, indexes, boundsType, bounds);

    return retVal;
}

