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

typedef struct
{
    ThreadSharedData *thShared;
    unsigned int rndState;
    int iterCount;

    int maxKickMagnitude;
    int *nodesToKick;

    Solution workingSol;

} ThreadSpecificData;

static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState)
{
    ThreadSpecificData thSpecific = {
        .thShared=thShared,
        .rndState=rndState,
        .iterCount=0
    };

    Instance *inst = thShared->bestSol.instance;
    thSpecific.workingSol=newSolution(inst);

    if (MAX_KICK_MAGNITUDE > inst->nNodes)
        thSpecific.maxKickMagnitude = inst->nNodes;

    // allocate more memory if permutations are necessary
    if (PERMUTATION_THRESHOLD < inst->nNodes)
        thSpecific.nodesToKick = malloc(inst->nNodes + thSpecific.maxKickMagnitude * sizeof(int));
    else
        thSpecific.nodesToKick = malloc((inst->nNodes - 1) * sizeof(int));
    
    if (!thSpecific.nodesToKick)
        throwError(inst, &thSpecific.workingSol, "NearestNeighbor -> initThreadSpecificData: Failed to allocate memory");

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    free(thSpecific->nodesToKick);
    thSpecific->nodesToKick = NULL;
}


static void *runVns(void* arg);

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

static void kick(ThreadSpecificData *thSpecific);

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

        if (!checkSolution(&thSpecific->workingSol))
            throwError(inst, &thSpecific->workingSol, "runVns: Solution after kick is incorrect (iteration = %d)", thSpecific->iterCount);

        apply2OptBestFix(&thSpecific->workingSol);

        if (thSpecific->workingSol.cost < thShared->bestSol.cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol.cost)
            {
                LOG(LOG_LVL_LOG, "Found better solution: cost = %f", thSpecific->workingSol.cost);
                Solution temp;
                swapElems(thSpecific->workingSol, thShared->bestSol, temp);
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

        int i = 0;
        while (i < currentKickMagnitude)
        {
            nodesToKick[i] = genRandom(&thSpecific->rndState, 1, n);

            // check if any is equal
            for (int j = 0; j < i; j++)
                if (nodesToKick[i] == nodesToKick[j])
                    continue;
            i++;
        }
    }

    // perform kick
    enum EdgeWeightType ewt = inst->params.edgeWeightType;
    bool roundW = inst->params.roundWeights;
    int *path = thSpecific->workingSol.indexPath;
    
    int first = path[nodesToKick[0]];

    for (int i = 0; i < currentKickMagnitude; i++)
    {
        // subtract old cost
        thSpecific->workingSol.cost -= (computeEdgeCost(inst->X[path[nodesToKick[i]-1]], inst->Y[path[nodesToKick[i]-1]], inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   ewt, roundW) +
                                        computeEdgeCost(inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   inst->X[path[nodesToKick[i]+1]], inst->Y[path[nodesToKick[i]+1]], ewt, roundW));

        int next;
        if (i < currentKickMagnitude-1)
            next = path[nodesToKick[i+1]];
        else
            next = first;

        path[nodesToKick[i]] = next;

        // add new cost
        thSpecific->workingSol.cost += computeEdgeCost(inst->X[path[nodesToKick[i]-1]], inst->Y[path[nodesToKick[i]-1]], inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   ewt, roundW) +
                                       computeEdgeCost(inst->X[path[nodesToKick[i]]],   inst->Y[path[nodesToKick[i]]],   inst->X[path[nodesToKick[i]+1]], inst->Y[path[nodesToKick[i]+1]], ewt, roundW);
    }
}

