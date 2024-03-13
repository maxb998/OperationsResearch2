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

static LocalBranchingData initLocalBranchingData (Solution *sol);

static void destroyLocalBranchingData (LocalBranchingData *lbData);

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

static void updateK(LocalBranchingData *lbData);

static void addSolSimilarityConstraint(LocalBranchingData *lbData);

static inline int removeSolSimilarityConstraint(LocalBranchingData *lbData);

void LocalBranching(Solution *sol, double timeLimit)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    LocalBranchingData lbData = initLocalBranchingData(sol);

    int errCode = CPXcallbacksetfunc(lbData.cpx.env, lbData.cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericCallbackCandidate, &lbData.cbData);
    if (errCode != 0)
        throwError("LocalBranching: CPXcallbacksetfunc failed with code %d", errCode);

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
        errCode = CPXsetdblparam(lbData.cpx.env, CPX_PARAM_TILIM, remainingTime);
        if (errCode != 0)
            throwError("LocalBranching: CPXsetdblparam failed with code %d", errCode);

        // update/set warm start solution
        LOG(LOG_LVL_TRACE, "LocalBranching: Setting WarmStart solution");
        errCode = WarmStart(&lbData.cpx, lbData.cbData.bestSuccessors);
        if (errCode != 0)
            throwError("LocalBranching: WarmStart failed with code %d", errCode);

        // run cplex
        LOG(LOG_LVL_NOTICE, "LocalBranching: Iteration %4d, running Branch&Cut method with k = %d", iterCount, lbData.k);
        LOG(LOG_LVL_TRACE, "LocalBranching: Running CPXmipopt");
        errCode = CPXmipopt(lbData.cpx.env, lbData.cpx.lp);
        if (errCode != 0)
            throwError("LocalBranching: CPXmipopt failed with code %d", errCode);

        errCode = removeSolSimilarityConstraint(&lbData);
        if (errCode != 0)
            throwError("LocalBranching: removeSolSimilarityConstraint -> CPXdelrows failed with code %d", errCode);

        if (lbData.k >= sol->instance->nNodes)
        {
            LOG(LOG_LVL_WARN, "LocalBranching: Solution similarity constraint is deactivated -> Optimal Solution found if finished before the time limit");
            break;
        }

        iterCount++;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }
    
    // update sol if necessary(very likely)
    if (lbData.cbData.bestCost < sol->cost)
        cvtSuccessorsToSolution(lbData.cbData.bestSuccessors, lbData.cbData.bestCost, sol);
    else
        LOG(LOG_LVL_WARN, "LocalBranching: Solution could not be optimized any further");

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    currentTime = cvtTimespec2Double(timeStruct);
    sol->execTime += currentTime - startTime;

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/sol->execTime);

    destroyLocalBranchingData(&lbData);
}

static void updateK(LocalBranchingData *lbData)
{
    static __uint128_t oldCost = -1LL;

    if (oldCost > lbData->cbData.bestCost)
    {
        oldCost = lbData->cbData.bestCost;
        lbData->k = K_MIN;
        return;
    }

    lbData->k += K_INCREMENT;
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
    if (errCode != 0)
        throwError("addSolSimilarityConstraint: CPXaddrows(): error %d", errCode);
}

static inline int removeSolSimilarityConstraint(LocalBranchingData *lbData)
{
    int n = lbData->sol->instance->nNodes;
    return CPXdelrows (lbData->cpx.env, lbData->cpx.lp, n, n);
}

