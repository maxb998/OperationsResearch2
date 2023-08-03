#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

#define MAX_KICK_MAGNITUDE 20
#define MIN_KICK_MAGNITUDE 5
#define PERMUTATION_THRESHOLD (MAX_KICK_MAGNITUDE * 2)
#define N_PERMUTATION(nNodes) (nNodes * 4)

typedef struct
{
    Solution bestSol;

    pthread_mutex_t mutex;

    double timeLimit;

} ThreadSharedData;

typedef struct
{
    ThreadSharedData *thShared;
    unsigned int rndState;
    int iterCount;

    int maxKickMagnitude;
    int *nodesToKick;

    Solution workingSol;
    float *X;
    float *Y;

} ThreadSpecificData;



static ThreadSharedData initThreadSharedData (Solution *sol, double timeLimit);
static void destroyThreadSharedData (ThreadSharedData *thShared);
static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState);
static void destroyThreadSpecificData(ThreadSpecificData *thSpecific);
static void *runVns(void* arg);
static void kick(ThreadSpecificData *thSpecific);


void VariableNeighborhoodSearch(Solution *sol, double timeLimit, int nThreads)
{
    Instance *inst = sol->instance;

    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // reset seed if debugging
    if (inst->params.logLevel == LOG_LVL_DEBUG)
        srand(inst->params.randomSeed);

    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError(sol->instance, sol, "VariableNeighborhood: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = inst->params.nThreads;

    if (!checkSolution(sol))
        throwError(inst, sol, "VariableNeighborhood: Input solution is not valid");

    sol->indexPath[inst->nNodes] = sol->indexPath[0];

    apply2OptBestFix(sol);

    ThreadSharedData thShared = initThreadSharedData(sol, startTime + timeLimit);
    ThreadSpecificData thSpecifics[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, (unsigned int)rand());
        pthread_create(&threads[i], NULL, runVns, &thSpecifics[i]);
    }

    int iterCount = 0;
    for (int i = 0; i < nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        iterCount += thSpecifics[i].iterCount;
        destroyThreadSpecificData(&thSpecifics[i]);
    }

    *sol = thShared.bestSol;

    destroyThreadSharedData(&thShared);    

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/sol->execTime);
}

static ThreadSharedData initThreadSharedData (Solution *sol, double timeLimit)
{
    ThreadSharedData thShared = { .timeLimit = timeLimit, .bestSol=*sol };

    if (pthread_mutex_init(&thShared.mutex, NULL)) throwError(sol->instance, sol, "VariableNeighborhoodSearch -> initThreadSharedData: Failed to initialize mutex");

    return thShared;
}

static void destroyThreadSharedData (ThreadSharedData *thShared)
{
    if (pthread_mutex_init(&thShared->mutex, NULL)) throwError(thShared->bestSol.instance, &thShared->bestSol, "VariableNeighborhoodSearch -> destroyThreadSharedData: Failed to destroy mutex");
}

static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState)
{
    ThreadSpecificData thSpecific = {
        .thShared=thShared,
        .rndState=rndState,
        .iterCount=0,
        .maxKickMagnitude=MAX_KICK_MAGNITUDE
    };

    Instance *inst = thShared->bestSol.instance;
    int n = inst->nNodes;

    thSpecific.workingSol=newSolution(inst);

    if (MAX_KICK_MAGNITUDE > n)
        thSpecific.maxKickMagnitude = n;

    // allocate more memory if permutations are necessary
    if (PERMUTATION_THRESHOLD < n)
        thSpecific.nodesToKick = malloc(thSpecific.maxKickMagnitude * sizeof(int));
    else
        thSpecific.nodesToKick = malloc((n - 1) * sizeof(int));
    
    if (!thSpecific.nodesToKick)
        throwError(inst, &thSpecific.workingSol, "NearestNeighbor -> initThreadSpecificData: Failed to allocate memory");
    
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        thSpecific.X = malloc((n + AVX_VEC_SIZE) * 2 * sizeof(int));
        if (!thSpecific.X)
            throwError(inst, &thSpecific.workingSol, "ExtraMileage -> initThreadSpecificData: Failed to allocate memory");
        thSpecific.Y = &thSpecific.X[n + AVX_VEC_SIZE];

        for (int i = 0; i < n; i++)
        {
            thSpecific.X[i] = inst->X[thShared->bestSol.indexPath[i]];
            thSpecific.Y[i] = inst->Y[thShared->bestSol.indexPath[i]];
        }
        thSpecific.X[n] = inst->X[thShared->bestSol.indexPath[0]];
        thSpecific.Y[n] = inst->Y[thShared->bestSol.indexPath[0]];
        for (int i = n+1; i < n + AVX_VEC_SIZE; i++)
        {
            thSpecific.X[i] = INFINITY;
            thSpecific.Y[i] = INFINITY;
        }
    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
        thSpecific.X = thSpecific.Y = NULL;
    #endif

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    free(thSpecific->nodesToKick);
    free(thSpecific->X);

    thSpecific->nodesToKick = NULL;
    thSpecific->X = NULL;
}

static void *runVns(void *arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    cloneSolution(&thShared->bestSol, &thSpecific->workingSol);

    // init nodesToKick if using permutations to kick
    if (inst->nNodes < PERMUTATION_THRESHOLD)
        for (int i = 0; i < inst->nNodes-1; i++)
            thSpecific->nodesToKick[i] = i+1;

    while (currentTime < thShared->timeLimit)
    {
        kick(thSpecific);

        LOG(LOG_LVL_EVERYTHING, "runVns: [%d] solution has been kicked. Cost=%f", thSpecific->workingSol.cost);

        if (!checkSolution(&thSpecific->workingSol))
            throwError(inst, &thSpecific->workingSol, "runVns: Solution after kick is incorrect (iteration = %d)", thSpecific->iterCount, thSpecific->iterCount);

        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->X, thSpecific->Y);
        #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
            apply2OptBestFix(&thSpecific->workingSol);
        #endif

        LOG(LOG_LVL_EVERYTHING, "runVns: [%d] solution has been optimized. Cost=%f", thSpecific->iterCount, thSpecific->workingSol.cost);

        if (thSpecific->workingSol.cost < thShared->bestSol.cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol.cost)
            {
                LOG(LOG_LVL_LOG, "Found better solution: cost = %lf", cvtCost2Double(thSpecific->workingSol.cost));
                cloneSolution(&thSpecific->workingSol, &thShared->bestSol);
            }
            pthread_mutex_unlock(&thShared->mutex);
        }

        thSpecific->iterCount++;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }
    
    return NULL;
}

static void kick(ThreadSpecificData *thSpecific)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    int *nodesToKick = thSpecific->nodesToKick;
    
    // decide how much this kick magnitude will be and which nodes it will affect(saved into thSpecific.nodesToKick)

    int currentKickMagnitude;

    if (n < PERMUTATION_THRESHOLD)
    {
        currentKickMagnitude = n-1;

        for (int i = 0; i < N_PERMUTATION(n); i++)
        {
            // get permutation indexes
            int rndIndex0 = genRandom(&thSpecific->rndState, 0, n-1), rndIndex1 = genRandom(&thSpecific->rndState, 0, n-1);
            while (rndIndex0 == rndIndex1) rndIndex1 = genRandom(&thSpecific->rndState, 0, n-1);

            register int temp;
            swapElems(nodesToKick[rndIndex0], nodesToKick[rndIndex1], temp);
        }
    }
    else
    {
        currentKickMagnitude = genRandom(&thSpecific->rndState, MIN_KICK_MAGNITUDE, MAX_KICK_MAGNITUDE);

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

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        enum EdgeWeightType ewt = inst->params.edgeWeightType;
        bool roundW = inst->params.roundWeights;
    #endif

    // must save somewhere the first value since it will be overwritten and won't be recovered otherwise
    int first = path[nodesToKick[0]];
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        float firstX = thSpecific->X[nodesToKick[0]];
        float firstY = thSpecific->Y[nodesToKick[0]];
    #endif
    

    for (int i = 0; i < currentKickMagnitude; i++)
    {
        float edgeCost0, edgeCost1;

        // subtract old cost
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            edgeCost0 = computeEdgeCost(thSpecific->X[nodesToKick[i]-1], thSpecific->Y[nodesToKick[i]-1], thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   ewt, roundW);
            edgeCost1 = computeEdgeCost(thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   thSpecific->X[nodesToKick[i]+1], thSpecific->Y[nodesToKick[i]+1], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            edgeCost0 = computeEdgeCost(inst->X[path[nodesToKick[i]-1]], inst->Y[path[nodesToKick[i]-1]], inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   ewt, roundW);
            edgeCost1 = computeEdgeCost(inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   inst->X[path[nodesToKick[i]+1]], inst->Y[path[nodesToKick[i]+1]], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            edgeCost0 = inst->edgeCostMat[path[nodesToKick[i]] * n + path[nodesToKick[i]-1]];
            edgeCost1 = inst->edgeCostMat[path[nodesToKick[i]] * n + path[nodesToKick[i]+1]];
        #endif

        thSpecific->workingSol.cost -= (cvtFloat2Cost(edgeCost0) + cvtFloat2Cost(edgeCost1));

        
        if (i < currentKickMagnitude-1)
        {
            path[nodesToKick[i]] = path[nodesToKick[i+1]];
            #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
                thSpecific->X[nodesToKick[i]] = thSpecific->X[nodesToKick[i+1]];
                thSpecific->Y[nodesToKick[i]] = thSpecific->Y[nodesToKick[i+1]];
            #endif
        }
        else
        {
            path[nodesToKick[i]] = first;
            #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
                thSpecific->X[nodesToKick[i]] = firstX;
                thSpecific->Y[nodesToKick[i]] = firstY;
            #endif
        }

        // add new cost
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            edgeCost0 = computeEdgeCost(thSpecific->X[nodesToKick[i]-1], thSpecific->Y[nodesToKick[i]-1], thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   ewt, roundW);
            edgeCost1 = computeEdgeCost(thSpecific->X[nodesToKick[i]],   thSpecific->Y[nodesToKick[i]],   thSpecific->X[nodesToKick[i]+1], thSpecific->Y[nodesToKick[i]+1], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            edgeCost0 = computeEdgeCost(inst->X[path[nodesToKick[i]-1]], inst->Y[path[nodesToKick[i]-1]], inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   ewt, roundW);
            edgeCost1 = computeEdgeCost(inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   inst->X[path[nodesToKick[i]+1]], inst->Y[path[nodesToKick[i]+1]], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            edgeCost0 = inst->edgeCostMat[path[nodesToKick[i]] * n + path[nodesToKick[i]-1]];
            edgeCost1 = inst->edgeCostMat[path[nodesToKick[i]] * n + path[nodesToKick[i]+1]];
        #endif
        
        thSpecific->workingSol.cost += (cvtFloat2Cost(edgeCost0) + cvtFloat2Cost(edgeCost1));
    }
}

