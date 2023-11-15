#include "TspCplex.h"
#include "Tsp.h"

#include <cplex.h>
#include <stdio.h>
#include <stdarg.h> // used for va_list
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// Amount of nodes of the instance below which no nodes will be fixed, effectively running branch&cut. Also in any point the algorithm fixAmount won't be going any lower than this
#define MIN_UNFIX 150
// Minimum amount of edges of the solution to fix(it doesn't make much sense fixing only 1 node, if nNodes = 151 then 50 nodes will be fixed and 101 will be free)
#define MIN_FIX 50
// Incremental/Decremental step in fixAmount during computation
#define FIX_OFFSET 10
// Number of non-improving iterations before fixAmount is increased
#define STATIC_COST_THRESHOLD 10LL

typedef struct
{
    CplexData cpx;

    CallbackData cbData;

    int fixAmount;

    // Bounds Value
    double *bVal;
    // Bounds Type
    char *bType;
    int *indexes;

} HardfixAllocatedMem;

// Initialize HardfixAllocatedMem struct allocating all necessary memory and initializing all sub-structs and initializes the fixAmount value
static HardfixAllocatedMem initHardfixAllocatedMem(Solution *sol);

// Free all memory allocated by init method and destriy all sub-structs
static void destroyHardfixAllocatedMem(HardfixAllocatedMem *hfAlloc);

// Updates fixAmount value. Only does so after STATIC_COST_THRESHOLD iterations
static void updateFixAmount(HardfixAllocatedMem *hfAlloc);

// Reset Upper and Lower bounds of the problem in cplex. Pass already allocated memory in indexes and bounds (ncols elements each)
static int resetBounds(HardfixAllocatedMem *hfAlloc);

// Fix edges in random way up the the fixAmount
static int randomFix(HardfixAllocatedMem *hfAlloc);

// Fix smallest edges up to the fixAmount
static int smallestFix(HardfixAllocatedMem *hfAlloc);


void HardFixing(Solution *sol, double timeLimit)
{
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTime = cvtTimespec2Double(currT);

    int errCode = 0;

    if (!checkSolution(sol))
		throwError("HardFixing: Input solution is not valid");

    HardfixAllocatedMem hfAlloc = initHardfixAllocatedMem(sol);

    if ((errCode = CPXcallbacksetfunc(hfAlloc.cpx.env, hfAlloc.cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericCallbackCandidate, &hfAlloc.cbData)) != 0)
        throwError("HardFix: CPXcallbacksetfunc failed with code %d", errCode);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double currentTime = cvtTimespec2Double(currT);

    int iterCount = 0;
    while (currentTime < startTime + timeLimit)
    {
        iterCount++;

        if (hfAlloc.fixAmount > 0)
        {
            switch (sol->instance->params.hardFixPolicy)
            {
            case HARDFIX_POLICY_RANDOM:
                if ((errCode = randomFix(&hfAlloc)) != 0)
                    throwError("HardFix: randomFix failed with code %d", errCode);
                break;
            case HARDFIX_POLICY_SMALLEST:
                if ((errCode = smallestFix(&hfAlloc)) != 0)
                    throwError("HardFix: smallestFix failed with code %d", errCode);
                break;
            }
        }

        #ifdef DEBUG
        LOG(LOG_LVL_ALL, "HardFix: Bounds were fixed");
        #endif

        // update time limit
        double remainingTime = startTime + timeLimit - currentTime;
        if ((errCode = CPXsetdblparam(hfAlloc.cpx.env, CPX_PARAM_TILIM, remainingTime)) != 0)
            throwError("HardFix: CPXsetdblparam failed with code %d", errCode);

        // update/set warm start solution
        if ((errCode = WarmStart(&hfAlloc.cpx, hfAlloc.cbData.bestSuccessors)) != 0)
            throwError("HardFix: WarmStart failed with code %d", errCode);

        #ifdef DEBUG
        LOG(LOG_LVL_ALL, "HardFix: WarmStart was set");
        #endif

        // run cplex
        LOG(LOG_LVL_DEBUG, "HardFix: Iteration %4d, running Branch&Cut method with %lu fixed edges", iterCount, hfAlloc.fixAmount);
        if ((errCode = CPXmipopt(hfAlloc.cpx.env, hfAlloc.cpx.lp)) != 0)
            throwError("HardFix: CPXmipopt failed with code %d", errCode);

        #ifdef DEBUG
        LOG(LOG_LVL_ALL, "HardFix: CPXmipopt finished");
        #endif

        if (hfAlloc.fixAmount == 0)
        {
            LOG(LOG_LVL_NOTICE, "HardFix: No edges were fixed. Best solution has been found");
            break;
        }
        
        updateFixAmount(&hfAlloc);

        if ((errCode = resetBounds(&hfAlloc)) != 0)
            throwError("HardFix: resetBounds failed with code %d", errCode);

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        currentTime = cvtTimespec2Double(currT);
    }
    
    // update sol if necessary(very likely)
    if (hfAlloc.cbData.bestCost < sol->cost)
        cvtSuccessorsToSolution(hfAlloc.cbData.bestSuccessors, sol);
    else
        LOG(LOG_LVL_WARNING, "HardFixing: Solution could not be optimized any further");

    destroyHardfixAllocatedMem(&hfAlloc);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    currentTime = cvtTimespec2Double(currT);

    sol->execTime += currentTime - startTime;

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/sol->execTime);
}

static HardfixAllocatedMem initHardfixAllocatedMem(Solution *sol)
{
    int n = sol->instance->nNodes;

    HardfixAllocatedMem hfAlloc = {
        .cpx = initCplexData(sol->instance),
        .bVal = malloc(n * sizeof(double)),
        .bType = malloc(n),
        .indexes = malloc(n * sizeof(int))
    };

    hfAlloc.cbData = initCallbackData(&hfAlloc.cpx, sol);

    // init fixAmount
    if(n < MIN_UNFIX)
    {
        LOG(LOG_LVL_WARNING, "HardFix: Solution is small, so no edge will be fixed, resulting in a computation that is the same as branch-cut");
        hfAlloc.fixAmount = 0;
    }
    else if (n < MIN_UNFIX + MIN_FIX)
        hfAlloc.fixAmount = MIN_FIX;
    else
        hfAlloc.fixAmount = n - MIN_UNFIX;

    return hfAlloc;
}

static void destroyHardfixAllocatedMem(HardfixAllocatedMem *hfAlloc)
{
    destroyCplexData(&hfAlloc->cpx);
    destroyCallbackData(&hfAlloc->cbData);
    free(hfAlloc->bVal);
    free(hfAlloc->bType);
    free(hfAlloc->indexes);

    hfAlloc->bVal = NULL;
    hfAlloc->bType = NULL;
    hfAlloc->indexes = NULL;
}

static void updateFixAmount(HardfixAllocatedMem *hfAlloc)
{
    // counts how for how many iterations the cost of the solution hasn't improved
    static int staticCostCount = -1;
    // keeps track of the previous iteration cost
    static __uint128_t oldCost = -1LL;

    if (oldCost <= hfAlloc->cbData.bestCost)
    {
        staticCostCount = 0;
        oldCost = hfAlloc->cbData.bestCost;
    }
    else
        staticCostCount++;
    
    if ((staticCostCount >= STATIC_COST_THRESHOLD) || ((staticCostCount > 0) && (hfAlloc->cpx.inst->params.hardFixPolicy == HARDFIX_POLICY_SMALLEST)))
    {
        hfAlloc->fixAmount -= FIX_OFFSET;

        if (hfAlloc->fixAmount < 0) // overflow of unsigned operation (fixAmount is way too big) -> set fixAmount to 0
            hfAlloc->fixAmount = 0;

        LOG(LOG_LVL_LOG, "Decreasing the amount of fixed edges to %u", hfAlloc->fixAmount);
        staticCostCount = 0; //(n - fixAmount - MIN_UNFIX) / FIX_OFFSET;
    }
    
}

static int resetBounds(HardfixAllocatedMem *hfAlloc)
{
    int n = hfAlloc->cpx.inst->nNodes;

    // set lower bounds (upper bounds are not modified)
    for (int i = 0; i < n; i++)
    {
        hfAlloc->bType[i] = 'L';
        hfAlloc->indexes[i] = xpos(i, hfAlloc->cbData.bestSuccessors[i], n);
        hfAlloc->bVal[i] = 0.0;
    }
    int retVal = CPXchgbds(hfAlloc->cpx.env, hfAlloc->cpx.lp, n, hfAlloc->indexes, hfAlloc->bType, hfAlloc->bVal);

    return retVal;
}

static int randomFix(HardfixAllocatedMem *hfAlloc)
{
    int n = hfAlloc->cpx.inst->nNodes;

    // setup indexes from 0 to n-1 in order to perform permutations on it
    for (int i = 0; i < n; i++)
        hfAlloc->indexes[i] = i;

    // n random permutations
    for (int i = 0; i < n; i++)
    {
        int pos1 = (int)((long)rand() * (long)hfAlloc->fixAmount / (long)RAND_MAX), pos2 = (int)((long)rand() * (long)hfAlloc->fixAmount / (long)RAND_MAX);
        while (pos2 == pos1) pos2 = (int)((long)rand() * (long)hfAlloc->fixAmount / (long)RAND_MAX);

        swapElems(hfAlloc->indexes[pos1], hfAlloc->indexes[pos2])
    }
    
    // setup arrays to fix bounds
    for (int i = 0; i < hfAlloc->fixAmount; i++)
    {
        hfAlloc->bType[i] = 'B';
        hfAlloc->indexes[i] = xpos(hfAlloc->indexes[i], hfAlloc->cbData.bestSuccessors[hfAlloc->indexes[i]], n);
        hfAlloc->bVal[i] = 1.0;
    }
    int retVal = CPXchgbds(hfAlloc->cpx.env, hfAlloc->cpx.lp, hfAlloc->fixAmount, hfAlloc->indexes, hfAlloc->bType, hfAlloc->bVal);

    return retVal;
}

static int smallestFix(HardfixAllocatedMem *hfAlloc)
{
    Instance *inst = hfAlloc->cpx.inst;
    int n = inst->nNodes;

    // instead of allocating new memory cast some pointers inside space allocated pointed by bounds(which is bounds)
    float *tourEdgesCost = (float*)hfAlloc->bVal;
    int *successors = hfAlloc->cbData.bestSuccessors;

    for (int i = n; i < n + AVX_VEC_SIZE; i++)
        successors[i] = 0;

    for (int i = 0; i < n; i++)
    {
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            tourEdgesCost[i] = computeEdgeCost(inst->X[i], inst->Y[i], inst->X[successors[i]], y2 = inst->Y[successors[i]], inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            tourEdgesCost[i] = inst->edgeCostMat[ i * n + successors[i]];
        #endif
    }
    
    # if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        // compute cost of each edge with avx instructions(not necessary since it's not a performance critical function)
        for (int i = 0; i < n; i += AVX_VEC_SIZE)
        {
            __m256 x1 = _mm256_loadu_ps(&inst->X[i]), y1 = _mm256_loadu_ps(&inst->Y[i]);

            __m256i indexesSucc = _mm256_loadu_si256((__m256i_u*)&successors[i]);
            __m256 x2 = _mm256_i32gather_ps(inst->X, indexesSucc, 4), y2 = _mm256_i32gather_ps(inst->Y, indexesSucc, 4);

            __m256 cost = computeEdgeCost_VEC(x1, y1, x2, y2, inst);

            _mm256_storeu_ps(&tourEdgesCost[i], cost);
        }
    # endif

    // perform argsort on tourEdgesCost using minCostIndexes as output
    argsort(tourEdgesCost, hfAlloc->indexes, n);

    for (int i = 0; i < hfAlloc->fixAmount; i++)
    {
        hfAlloc->bType[i] = 'B';
        hfAlloc->indexes[i] = xpos(hfAlloc->indexes[i], hfAlloc->cbData.bestSuccessors[hfAlloc->indexes[i]], n);
        hfAlloc->bVal[i] = 1.0;
    }
    
    int retVal = CPXchgbds(hfAlloc->cpx.env, hfAlloc->cpx.lp, hfAlloc->fixAmount, hfAlloc->indexes, hfAlloc->bType, hfAlloc->bVal);

    return retVal;
}

