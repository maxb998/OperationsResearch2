#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

#define PERMUTATION_THRESHOLD(maxKick) (maxKick * 2)
#define N_PERMUTATION(nNodes) (nNodes * 4)

typedef struct
{
    Solution *bestSol;

    pthread_mutex_t mutex;

    double timeLimit;

} ThreadSharedData;

typedef struct
{
    ThreadSharedData *thShared;
    unsigned int rndState;
    int iterCount;

    // pointer to memory allocation where selected nodes to be kicked are stored
    int *nodesToKick;

    Solution workingSol;
    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
    float *X;
    float *Y;
    #endif
    float *costCache;

} ThreadSpecificData;


// Initializes mutex and variables
static ThreadSharedData initThreadSharedData (Solution *sol, double timeLimit);
// Destroy mutex
static void destroyThreadSharedData (ThreadSharedData *thShared);
// Initializes variables (and allocates memory on .X, .Y when using COMPUTE_OPTION_AVX)
static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState);
// Frees memory allocated in X and Y if any
static void destroyThreadSpecificData(ThreadSpecificData *thSpecific);
// Method each thread run when executing VNS until time limit
static void *runVns(void* arg);
// Set the workingSol equal to thShared.bestSol (and rearranges .X, .Y accordingly when using COMPUTE_OPTION_AVX)
static void setupThSpecificOnBestSol(ThreadSpecificData *thSpecific);
// Give a random Kick to the solution by swapping points inside thSpecific.workingSol (and .X, .Y when using COMPUTE_OPTION_AVX)
static void kick(ThreadSpecificData *thSpecific);



void VariableNeighborhoodSearch(Solution *sol, double timeLimit)
{
    Instance *inst = sol->instance;

    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // must save solution time since 2opt is gonna increase it inconsistently with vns
    double initialSolutionRuntime = sol->execTime;

    // reset seed if debugging
    if (inst->params.logLevel == LOG_LVL_DEBUG)
        srand(inst->params.randomSeed);

    if (!checkSolution(sol))
        throwError("VariableNeighborhood: Input solution is not valid");

    sol->indexPath[inst->nNodes] = sol->indexPath[0];

    apply2OptBestFix(sol); // isn't necessary since it's solution should already be 2-optimized, but just to be sure

    ThreadSharedData thShared = initThreadSharedData(sol, startTime + timeLimit);
    ThreadSpecificData thSpecifics[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, (unsigned int)rand());
        pthread_create(&threads[i], NULL, runVns, &thSpecifics[i]);
    }

    int iterCount = 0;
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        iterCount += thSpecifics[i].iterCount;
        destroyThreadSpecificData(&thSpecifics[i]);
    }

    destroyThreadSharedData(&thShared);   

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double execTime = cvtTimespec2Double(timeStruct) - startTime;
    sol->execTime = execTime + initialSolutionRuntime;

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %d", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/execTime);
}

static ThreadSharedData initThreadSharedData (Solution *sol, double timeLimit)
{
    ThreadSharedData thShared = { .timeLimit = timeLimit, .bestSol=sol };

    if (pthread_mutex_init(&thShared.mutex, NULL)) throwError("VariableNeighborhoodSearch -> initThreadSharedData: Failed to initialize mutex");

    return thShared;
}

static void destroyThreadSharedData (ThreadSharedData *thShared)
{
    if (pthread_mutex_init(&thShared->mutex, NULL)) throwError("VariableNeighborhoodSearch -> destroyThreadSharedData: Failed to destroy mutex");
}

static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState)
{
    ThreadSpecificData thSpecific = {
        .thShared=thShared,
        .rndState=rndState,
        .iterCount=0
    };

    Instance *inst = thShared->bestSol->instance;
    int n = inst->nNodes;

    thSpecific.workingSol=newSolution(inst);

    if (inst->params.vnsKickSize.Max > n)
        throwError("Vns maximum kick cannot be greater than the size of the instance, at most equal");

    size_t memToAlloc_kick;
    if (PERMUTATION_THRESHOLD(inst->params.vnsKickSize.Max) < n)
        memToAlloc_kick = inst->params.vnsKickSize.Max * sizeof(int);
    else
        memToAlloc_kick = (n - 1) * sizeof(int);

    size_t memToAlloc_other = n + AVX_VEC_SIZE;
    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        memToAlloc_other += (n + AVX_VEC_SIZE) * 2;
    #endif
    memToAlloc_other *= sizeof(float);

    thSpecific.nodesToKick = malloc(memToAlloc_kick + memToAlloc_other);
    if (!thSpecific.nodesToKick)
        throwError("VariableNeighborhoodSearch -> initThreadSpecificData: Failed to allocate memory");
    thSpecific.costCache = (float*)&thSpecific.nodesToKick[memToAlloc_kick / sizeof(int)];
    
    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        thSpecific.X = &thSpecific.costCache[n + AVX_VEC_SIZE];
        thSpecific.Y = &thSpecific.X[n + AVX_VEC_SIZE];
    #endif

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    free(thSpecific->nodesToKick);

    thSpecific->nodesToKick = NULL;
    thSpecific->costCache = NULL;
    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        thSpecific->X = thSpecific->Y = NULL;
    #endif
}

static void *runVns(void *arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    setupThSpecificOnBestSol(thSpecific);

    // init nodesToKick if using permutations to kick
    if (inst->nNodes < PERMUTATION_THRESHOLD(inst->params.vnsKickSize.Max))
        for (int i = 0; i < inst->nNodes-1; i++)
            thSpecific->nodesToKick[i] = i+1;

    int nonImprovingIterCount = 0;
    int restartThreshold = thSpecific->workingSol.instance->params.metaRestartThreshold;

    while (currentTime < thShared->timeLimit)
    {
        // after some non-improving iterations set the workingSol equal to the bestSol
        if (nonImprovingIterCount > restartThreshold)
        {
            setupThSpecificOnBestSol(thSpecific);
            nonImprovingIterCount = 0;
        }

        kick(thSpecific);

        LOG(LOG_LVL_EVERYTHING, "runVns: [%d] solution has been kicked. Cost=%lf", cvtCost2Double(thSpecific->workingSol.cost));

        #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->X, thSpecific->Y, thSpecific->costCache);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->costCache);
        #endif

        LOG(LOG_LVL_EVERYTHING, "runVns: [%d] solution has been optimized. Cost=%lf", thSpecific->iterCount, cvtCost2Double(thSpecific->workingSol.cost));

        if (thSpecific->workingSol.cost < thShared->bestSol->cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol->cost)
            {
                LOG(LOG_LVL_LOG, "Found better solution: cost = %lf", cvtCost2Double(thSpecific->workingSol.cost));
                cloneSolution(&thSpecific->workingSol, thShared->bestSol);
            }
            else
                nonImprovingIterCount++;
            pthread_mutex_unlock(&thShared->mutex);
        }
        else
            nonImprovingIterCount++;

        thSpecific->iterCount++;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }
    
    return NULL;
}

static void setupThSpecificOnBestSol(ThreadSpecificData *thSpecific)
{
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    // clone of best sol must be done with locked mutex otherwise there may be synchronization errors when starting a lot of threads
    pthread_mutex_lock(&thShared->mutex);
    cloneSolution(thShared->bestSol, &thSpecific->workingSol);
    pthread_mutex_unlock(&thShared->mutex);

    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        for (int i = 0; i < n; i++)
        {
            thSpecific->X[i] = inst->X[thSpecific->workingSol.indexPath[i]];
            thSpecific->Y[i] = inst->Y[thSpecific->workingSol.indexPath[i]];
        }
        thSpecific->X[n] = thSpecific->X[0];
        thSpecific->Y[n] = thSpecific->Y[0];

        for (int i = n+1; i < n + AVX_VEC_SIZE; i++)
        {
            thSpecific->X[i] = INFINITY;
            thSpecific->Y[i] = INFINITY;
        }
    #endif

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        for (int i = 0; i < n; i++) // build cost cache
            thSpecific->costCache[i] = computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[i+1], thSpecific->Y[i+1], inst);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        int *path = thShared->bestSol->indexPath;
        for (int i = 0; i < n; i++) // build cost cache
            thSpecific->costCache[i] = inst->edgeCostMat[path[i] * n + path[i+1]];
    #endif

    for (int i = n+1; i < n + AVX_VEC_SIZE; i++)
        thSpecific->costCache[i] = INFINITY;

}

static void kick(ThreadSpecificData *thSpecific)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    int *nodesToKick = thSpecific->nodesToKick;
    
    // decide how much this kick magnitude will be and which nodes it will affect(saved into thSpecific.nodesToKick)

    int currentKickMagnitude;

    if (n < PERMUTATION_THRESHOLD(inst->params.vnsKickSize.Max))
    {
        currentKickMagnitude = n-1;

        for (int i = 0; i < N_PERMUTATION(n); i++)
        {
            // get permutation indexes
            int rndIndex0 = genRandom(&thSpecific->rndState, 0, n-1), rndIndex1 = genRandom(&thSpecific->rndState, 0, n-1);
            while (rndIndex0 == rndIndex1) rndIndex1 = genRandom(&thSpecific->rndState, 1, n-1);

            swapElems(nodesToKick[rndIndex0], nodesToKick[rndIndex1])
        }
    }
    else
    {
        currentKickMagnitude = genRandom(&thSpecific->rndState, inst->params.vnsKickSize.Min, inst->params.vnsKickSize.Max+1);

        for (int i = 0; i < currentKickMagnitude; i++)
        {
            nodesToKick[i] = genRandom(&thSpecific->rndState, 1, n);

            // check if any is equal
            for (int j = 0; j < i; j++)
            {
                if (nodesToKick[i] == nodesToKick[j])
                {
                    i--;
                    break;
                }
            }
        }
    }

    // perform kick
    int *path = thSpecific->workingSol.indexPath;

    // must save somewhere the first value since it will be overwritten and won't be recovered otherwise
    int first = path[nodesToKick[0]];
    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        float firstX = thSpecific->X[nodesToKick[0]];
        float firstY = thSpecific->Y[nodesToKick[0]];
    #endif
    

    for (int i = 0; i < currentKickMagnitude; i++)
    {
        float edgeCost0, edgeCost1;

        // subtract old cost
        #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            edgeCost0 = computeEdgeCost(thSpecific->X[nodesToKick[i]-1], thSpecific->Y[nodesToKick[i]-1], thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   inst);
            edgeCost1 = computeEdgeCost(thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   thSpecific->X[nodesToKick[i]+1], thSpecific->Y[nodesToKick[i]+1], inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            edgeCost0 = inst->edgeCostMat[(size_t)path[nodesToKick[i]] * (size_t)n + (size_t)path[nodesToKick[i]-1]];
            edgeCost1 = inst->edgeCostMat[(size_t)path[nodesToKick[i]] * (size_t)n + (size_t)path[nodesToKick[i]+1]];
        #endif

        thSpecific->workingSol.cost -= (cvtFloat2Cost(edgeCost0) + cvtFloat2Cost(edgeCost1));

        
        if (i < currentKickMagnitude-1)
        {
            path[nodesToKick[i]] = path[nodesToKick[i+1]];
            #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                thSpecific->X[nodesToKick[i]] = thSpecific->X[nodesToKick[i+1]];
                thSpecific->Y[nodesToKick[i]] = thSpecific->Y[nodesToKick[i+1]];
            #endif
        }
        else
        {
            path[nodesToKick[i]] = first;
            #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                thSpecific->X[nodesToKick[i]] = firstX;
                thSpecific->Y[nodesToKick[i]] = firstY;
            #endif
        }

        // add new cost
        #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            edgeCost0 = computeEdgeCost(thSpecific->X[nodesToKick[i]-1], thSpecific->Y[nodesToKick[i]-1], thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   inst);
            edgeCost1 = computeEdgeCost(thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   thSpecific->X[nodesToKick[i]+1], thSpecific->Y[nodesToKick[i]+1], inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            edgeCost0 = inst->edgeCostMat[(size_t)path[nodesToKick[i]] * (size_t)n + (size_t)path[nodesToKick[i]-1]];
            edgeCost1 = inst->edgeCostMat[(size_t)path[nodesToKick[i]] * (size_t)n + (size_t)path[nodesToKick[i]+1]];
        #endif
        
        thSpecific->workingSol.cost += (cvtFloat2Cost(edgeCost0) + cvtFloat2Cost(edgeCost1));
        thSpecific->costCache[nodesToKick[i]-1] = edgeCost0;
        thSpecific->costCache[nodesToKick[i]] = edgeCost1;
    }

    if ((inst->params.logLevel >= LOG_LVL_DEBUG) && (!checkSolution(&thSpecific->workingSol)))
    {
        for (int i = 0; i < currentKickMagnitude; i++)
            LOG(LOG_LVL_DEBUG, "iter n° %d nodesToKick[%d] = %d", thSpecific->iterCount, i, nodesToKick[i]);
        throwError("runVns: Solution after kick is incorrect (iteration = %d)", thSpecific->iterCount);
    }
}

