#include "TspCplex.h"
#include "Tsp.h"

#include <cplex.h>
#include <stdio.h>
#include <stdarg.h> // used for va_list
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


#define START_FIX_AMOUNT 0.9
// Percentage of
#define OFFSET_FIX_AMOUNT 0.1
// Number of non-improving iterations before fixAmount is increased
#define STATIC_COST_THRESHOLD 10

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

// Fix a random consecutive sequence of lenght fixAmount
static int sequenceFix(HardfixAllocatedMem *hfAlloc);


void HardFixing(Solution *sol, double timeLimit)
{
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTime = cvtTimespec2Double(currT);

    int errCode = 0;

    if (!checkSolution(sol))
		throwError("HardFixing: Input solution is not valid");

    HardfixAllocatedMem hfAlloc = initHardfixAllocatedMem(sol);

    errCode = CPXcallbacksetfunc(hfAlloc.cpx.env, hfAlloc.cpx.lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericCallbackCandidate, &hfAlloc.cbData);
    if (errCode != 0)
        throwError("HardFix: CPXcallbacksetfunc failed with code %d", errCode);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double currentTime = cvtTimespec2Double(currT);

    int iterCount = 0;
    while (currentTime < startTime + timeLimit)
    {
        if (hfAlloc.fixAmount > 0)
        {
            int rndNum = (long)rand() * 3 / RAND_MAX;
            switch (rndNum)
            {
            case 0:
                LOG(LOG_LVL_DEBUG, "Using randomFix at iteration %d", iterCount);
                errCode = randomFix(&hfAlloc);
                if (errCode != 0)
                    throwError("HardFix: randomFix failed with code %d", errCode);
                break;
            case 1:
                LOG(LOG_LVL_DEBUG, "Using smallestFix at iteration %d", iterCount);
                errCode = smallestFix(&hfAlloc);
                if (errCode != 0)
                    throwError("HardFix: smallestFix failed with code %d", errCode);
                break;
            default: // or case 2
                LOG(LOG_LVL_DEBUG, "Using sequenceFix at iteration %d", iterCount);
                errCode = sequenceFix(&hfAlloc);
                if (errCode != 0)
                    throwError("HardFix: sequenceFix failed with code %d", errCode);
                break;
            }
        }

        #ifdef DEBUG
        LOG(LOG_LVL_TRACE, "HardFix: Bounds were fixed");
        #endif

        // update time limit
        double remainingTime = startTime + timeLimit - currentTime;
        errCode = CPXsetdblparam(hfAlloc.cpx.env, CPX_PARAM_TILIM, remainingTime);
        if (errCode != 0)
            throwError("HardFix: CPXsetdblparam failed with code %d", errCode);

        // update/set warm start solution
        errCode = WarmStart(&hfAlloc.cpx, hfAlloc.cbData.bestSuccessors);
        if (errCode != 0)
            throwError("HardFix: WarmStart failed with code %d", errCode);

        #ifdef DEBUG
        LOG(LOG_LVL_TRACE, "HardFix: WarmStart was set");
        #endif

        // run cplex
        LOG(LOG_LVL_DEBUG, "HardFix: Iteration %4d, running Branch&Cut method with %lu fixed edges", iterCount, hfAlloc.fixAmount);
        errCode = CPXmipopt(hfAlloc.cpx.env, hfAlloc.cpx.lp);
        if (errCode != 0)
            throwError("HardFix: CPXmipopt failed with code %d", errCode);

        #ifdef DEBUG
        LOG(LOG_LVL_TRACE, "HardFix: CPXmipopt finished");
        #endif

        if (hfAlloc.fixAmount == 0)
        {
            int solStatus = CPXgetstat(hfAlloc.cpx.env, hfAlloc.cpx.lp);
            if (((solStatus == CPXMIP_OPTIMAL) || (solStatus == CPXMIP_OPTIMAL_TOL)))
                LOG(LOG_LVL_NOTICE, "Optimal solution found");
            break;
        }
        
        updateFixAmount(&hfAlloc);

        errCode = resetBounds(&hfAlloc);
        if (errCode != 0)
            throwError("HardFix: resetBounds failed with code %d", errCode);

        iterCount++;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        currentTime = cvtTimespec2Double(currT);
    }
    
    // update sol if necessary(very likely)
    if (hfAlloc.cbData.bestCost < sol->cost)
        cvtSuccessorsToSolution(hfAlloc.cbData.bestSuccessors, hfAlloc.cbData.bestCost, sol);
    else
        LOG(LOG_LVL_WARN, "HardFixing: Solution could not be optimized any further");

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
        .bVal = malloc(n * (sizeof(double) + sizeof(int) + sizeof(char)))
    };
    hfAlloc.bType = (char*)&hfAlloc.bVal[n];
    hfAlloc.indexes = (int*)&hfAlloc.bType[n];

    hfAlloc.cbData = initCallbackData(&hfAlloc.cpx, sol);

    // init fixAmount
    hfAlloc.fixAmount = (int)((float)n * START_FIX_AMOUNT);

    return hfAlloc;
}

static void destroyHardfixAllocatedMem(HardfixAllocatedMem *hfAlloc)
{
    destroyCplexData(&hfAlloc->cpx);
    destroyCallbackData(&hfAlloc->cbData);
    free(hfAlloc->bVal);

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

    if (oldCost > hfAlloc->cbData.bestCost)
    {
        staticCostCount = 0;
        oldCost = hfAlloc->cbData.bestCost;
    }
    else
        staticCostCount++;
    
    if (staticCostCount >= STATIC_COST_THRESHOLD)
    {
        hfAlloc->fixAmount = hfAlloc->fixAmount - (int)ceilf(((float)hfAlloc->cbData.inst->nNodes * OFFSET_FIX_AMOUNT));

        if (hfAlloc->fixAmount < 0) // overflow of unsigned operation (fixAmount is way too big) -> set fixAmount to 0
            hfAlloc->fixAmount = 0;

        LOG(LOG_LVL_INFO, "Decreasing the amount of fixed edges to %u", hfAlloc->fixAmount);
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
        int pos1 = (int)((long)rand() * n / RAND_MAX), pos2 = (int)((long)rand() * n / RAND_MAX);
        while (pos2 == pos1) pos2 = (int)((long)rand() * n / RAND_MAX);

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

    for (int i = n; i < 2*n + AVX_VEC_SIZE; i++)
        successors[i] = 0;

    for (int i = 0; i < n; i++)
    {
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            tourEdgesCost[i] = computeEdgeCost(inst->X[i], inst->Y[i], inst->X[successors[i]], inst->Y[successors[i]], inst);
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

static int sequenceFix(HardfixAllocatedMem *hfAlloc)
{
    int n = hfAlloc->cpx.inst->nNodes;

    int seqStartPt = (int)((long)rand() * n / RAND_MAX);
    
    // setup arrays to fix bounds
    for (int i = 0, j = seqStartPt; i < hfAlloc->fixAmount; i++, j=hfAlloc->cbData.bestSuccessors[j])
    {   
        hfAlloc->bType[i] = 'B';
        hfAlloc->indexes[i] = xpos(j, hfAlloc->cbData.bestSuccessors[j], n);
        hfAlloc->bVal[i] = 1.0;
    }
    int retVal = CPXchgbds(hfAlloc->cpx.env, hfAlloc->cpx.lp, hfAlloc->fixAmount, hfAlloc->indexes, hfAlloc->bType, hfAlloc->bVal);

    return retVal;
}
