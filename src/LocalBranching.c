#include "TspCplex.h"
#include "Tsp.h"

#include <cplex.h>
#include <stdio.h>
#include <stdarg.h> // used for va_list
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

#define K_MIN 10
#define K_INCREMENT 5

typedef struct 
{
    CplexData cpx;
    CallbackData cbData;
    Solution *sol;
    int k;
} LocalBranchingData;

static LocalBranchingData initLocalBranchingData (Solution *sol)
{
    LocalBranchingData lbData = {
        .cpx = initCplexData(sol->instance),
        .cbData = initCallbackData(&lbData.cpx, sol),
        .sol = sol,
        .k = K_MIN
    };

    return lbData;
}

static void destroyLocalBranchingData (LocalBranchingData *lbData)
{
    destroyCallbackData(&lbData->cbData);
    destroyCplexData(&lbData->cpx);
}

static void throwLocalBranchingError(LocalBranchingData *lbData, Solution *sol, char *msg, ...)
{
    printf("[\033[1;31mERR \033[0m] ");

    va_list params;
    va_start(params, msg);
    vprintf(msg, params);
    va_end(params);

    printf("\n");

    if (lbData) destroyLocalBranchingData(lbData);
    if (sol->instance) destroyInstance(sol->instance);
    if (sol) destroySolution(sol);

    exit(EXIT_FAILURE);
}


static void updateK(LocalBranchingData *lbData);

static void addSolSimilarityConstraint(LocalBranchingData *lbData);

static inline int removeSolSimilarityConstraint(LocalBranchingData *lbData);

void LocalBranching(Solution *sol, double timeLimit)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    LocalBranchingData lbData = initLocalBranchingData(sol);

    int errCode;

    if ((errCode = CPXcallbacksetfunc(lbData.cpx.env, lbData.cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericCallbackCandidate, &lbData.cbData)))
        throwLocalBranchingError(&lbData, sol, "LocalBranching: CPXcallbacksetfunc failed with code %d", errCode);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    int iterCount = 0;
    while (currentTime < startTime + timeLimit)
    {
        updateK(&lbData);

        // add constraint
        addSolSimilarityConstraint(&lbData);

        // update time limit
        double remainingTime = startTime + timeLimit - currentTime;
        if ((errCode = CPXsetdblparam(lbData.cpx.env, CPX_PARAM_TILIM, remainingTime)))
            throwLocalBranchingError(&lbData, sol, "LocalBranching: CPXsetdblparam failed with code %d", errCode);

        // update/set warm start solution
        if ((errCode = WarmStart(&lbData.cpx, lbData.cbData.bestSuccessors)))
            throwLocalBranchingError(&lbData, sol, "LocalBranching: WarmStart failed with code %d", errCode);

        LOG(LOG_LVL_DEBUG, "LocalBranching: WarmStart was set");

        // run cplex
        LOG(LOG_LVL_NOTICE, "LocalBranching: Iteration %4d, running Branch&Cut method with k = %d", iterCount, lbData.k);
        if ((errCode = CPXmipopt(lbData.cpx.env, lbData.cpx.lp)))
            throwLocalBranchingError(&lbData, sol, "LocalBranching: CPXmipopt failed with code %d", errCode);

        LOG(LOG_LVL_DEBUG, "LocalBranching: CPXmipopt finished");

        if ((errCode = removeSolSimilarityConstraint(&lbData)))
            throwLocalBranchingError(&lbData, sol, "LocalBranching: removeSolSimilarityConstraint -> CPXdelrows failed with code %d", errCode);

        if (lbData.k >= sol->instance->nNodes)
        {
            LOG(LOG_LVL_WARNING, "LocalBranching: Solution similarity constraint is deactivated -> Optimal Solution found");
            break;
        }

        iterCount++;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }
    
    // update sol if necessary(very likely)
    if (lbData.cbData.bestCost < sol->cost)
    {
        cvtSuccessorsToSolution(lbData.cbData.bestSuccessors, sol);
        sol->cost = computeSolutionCost(sol);
    }
    else
        LOG(LOG_LVL_WARNING, "LocalBranching: Solution could not be optimized any further");

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    currentTime = cvtTimespec2Double(timeStruct);
    sol->execTime += currentTime - startTime;

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/sol->execTime);

    destroyLocalBranchingData(&lbData);
}

static void updateK(LocalBranchingData *lbData)
{
    static double oldCost = INFINITY;

    if (oldCost > lbData->cbData.bestCost)
    {
        oldCost = lbData->cbData.bestCost;
        lbData->k = K_MIN;
        return;
    }

    lbData->k += K_INCREMENT;

    //if (lbData->k)
}

static void addSolSimilarityConstraint(LocalBranchingData *lbData)
{
    int n = lbData->sol->instance->nNodes;
    int *successors = lbData->cbData.bestSuccessors;

    int *index = malloc(n * sizeof(int));
    double *value = malloc(n * sizeof(double));

    double rhs = n - lbData->k;

    for (int i = 0; i < n; i++)
    {
        index[i] = xpos(i, successors[i], n);
        value[i] = 1.0;
    }

    const char sense = 'G'; // 'G' for Greater or Equal
    char *rname = "Local Branching constraint";
    const int izero = 0;

    int errCode = CPXaddrows(lbData->cpx.env, lbData->cpx.lp, 0, 1, n, &rhs, &sense, &izero, index, value, NULL, &rname);
    if (errCode)
    {
        destroyLocalBranchingData(lbData); free(index); free(value);
        throwError(lbData->sol->instance, lbData->sol, "addSolSimilarityConstraint: CPXaddrows(): error %d", errCode);
    }
}

static inline int removeSolSimilarityConstraint(LocalBranchingData *lbData)
{
    int n = lbData->sol->instance->nNodes;
    return CPXdelrows (lbData->cpx.env, lbData->cpx.lp, n, n);
}

