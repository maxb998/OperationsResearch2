#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

typedef struct
{
    double offset;
    int index1;
    int index2;

}SwapInformation;

typedef struct
{
    Solution *bestSol;
    pthread_mutex_t mutex;
    double startingTemperature;
    double startTime;
    double timeLimit;

} ThreadSharedData;

typedef struct
{
    ThreadSharedData *thShared;
    Solution thSol;
    unsigned int rndState;
    int iterations;
    int iterationsSinceUpdate;
    double threshold;
    double temperature;

} ThreadSpecificData;

static inline ThreadSharedData initThreadSharedData(Solution *sol, double timeLimit, double startTime);
static inline void destroyThreadSharedData(ThreadSharedData *thShared);
static inline ThreadSpecificData initThreadSpecificData(Solution *sol, ThreadSharedData *thShared, unsigned int rndState);
static inline void destroyThreadSpecificData(ThreadSpecificData *thSpecific);

static inline SwapInformation randomSwap(Solution *sol, unsigned int * rndState);

static inline void updateTemperature(float *temperature);

static inline void normalizeDelta(double *deltaCost, Solution *sol);

static bool keepMove(double threshold);

void SimulatedAnnealing(Solution *sol, double timeLimit)
{

    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // check of the input solution
    if (!checkSolution(sol))
    throwError(sol->instance, sol, "SimulatedAnnealing: Input solution is not valid");

    // initialization of threads data
    ThreadSharedData thShared = initThreadSharedData(sol, timeLimit, startTime);

    ThreadSpecificData thSpecific[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < sol->instance->params.nThreads; i++)
    {
        thSpecific[i] = initThreadSpecificData(sol, &thShared, (unsigned int)rand());
        pthread_create(&threads[i], NULL, runSimulatedAnnealing, &thSpecific[i]);
    }

    int iterations = 0;
    for (int i = 0; i < sol->instance->params.nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        iterations += thSpecific[i].iterations;
    }

    cloneSolution(thShared.bestSol, &sol);
    checkSolution(&sol);

    for (int i = 0; i < sol->instance->params.nThreads; i++)
        destroyThreadSpecificData(&thSpecific[i]);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);
    LOG(LOG_LVL_NOTICE, "Total number of iterations: %d", iterations);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterations/(currentTime-startTime));
}

static inline ThreadSharedData initThreadSharedData(Solution *sol, double timeLimit, double startTime)
{
    ThreadSharedData thShared = { 
        .bestSol = sol,
        .timeLimit = startTime + timeLimit,
        .startTime = startTime
    };
    if(sol->instance->params.annealingTemperature > 0) thShared.startingTemperature = sol->instance->params.annealingTemperature;
    else thShared.startingTemperature = 10000;
    pthread_mutex_init(&thShared.mutex, NULL);
    return thShared;
}

static inline void destroyThreadSharedData(ThreadSharedData *thShared)
{
    pthread_mutex_destroy(&thShared->mutex);
}

static inline ThreadSpecificData initThreadSpecificData(Solution *sol, ThreadSharedData *thShared, unsigned int rndState)
{
    ThreadSpecificData thSpecific = {
        .iterations = 0,
        .iterationsSinceUpdate = 0,
        .rndState = rndState,
        .threshold = 0.0,
        .thShared = thShared,
        .thSol = NULL,
        .temperature = thShared->startingTemperature,
    };
    cloneSolution(sol, &thSpecific.thSol);

    return thSpecific;
}

static inline void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{

}

static void * runSimulatedAnnealing(void * arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;

    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);
 
    SwapInformation swapInfo = {
        .index1 = -1,
        .index2 = -1,
        .offset = 0
    };  

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    while (currentTime < thShared->timeLimit)
    {
        thSpecific->iterations++;
        thSpecific->iterationsSinceUpdate++;
        swapInfo = randomSwap(&thSpecific->thSol, thSpecific->rndState);

        if (swapInfo.offset < 0) // then we keep the swap
        {
            swapElems(thSpecific->thSol.indexPath[swapInfo.index1], thSpecific->thSol.indexPath[swapInfo.index2]);
            thSpecific->thSol.cost = computeSolutionCost(&thSpecific->thSol);
            LOG(LOG_LVL_LOG, "Better solution found at iteration: %d\t new cost: %lf", thSpecific->iterations, cvtCost2Double(thSpecific->thSol.cost));
            thSpecific->iterationsSinceUpdate = 0;
        }
        else // we implement the move with some probability
        {
            normalizeDelta(&swapInfo.offset, &thSpecific->thSol);
            thSpecific->iterations = exp(-swapInfo.offset/thSpecific->temperature);

            if (keepMove(thSpecific->threshold)) // then we keep the swap
            {
                swapElems(thSpecific->thSol.indexPath[swapInfo.index1], thSpecific->thSol.indexPath[swapInfo.index2]);
                thSpecific->thSol.cost = computeSolutionCost(&thSpecific->thSol);
                LOG(LOG_LVL_EVERYTHING, "Accepting bad move at iteration: %d\t, new cost: %lf", thSpecific->iterations, cvtCost2Double(thSpecific->thSol.cost));
                thSpecific->iterationsSinceUpdate = 0;
            }
            
            updateTemperature(&thSpecific->temperature);
        }

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);

        // if newSol has better cost than sol we update sol
        pthread_mutex_lock(&thShared->mutex);
        if(thSpecific->thSol.cost < thShared->bestSol->cost) 
            cloneSolution(&thSpecific->thSol, thShared->bestSol);
            thShared->bestSol->execTime += currentTime - thShared->startTime;
        pthread_mutex_unlock(&thShared->mutex);


        // if solution has not been updated for metaRestartThreshold iterations, and the time limit hasn't passed yet, we restart the annealing process from the best solution found so far
        if(thSpecific->iterationsSinceUpdate == thShared->bestSol->instance->params.metaRestartThreshold)
        {
            LOG(LOG_LVL_DEBUG, "Restarting Simulated Annealing from best solution found so far");
            thSpecific->iterationsSinceUpdate = 0;
            thSpecific->temperature = thShared->startingTemperature;
            cloneSolution(thShared->bestSol, &thSpecific->thSol);

        }
        
    }
}

static inline SwapInformation randomSwap(Solution *sol, unsigned int * rndState)
{
    int nNodes = sol->instance->nNodes;
    int *X = sol->instance->X;
    int *Y = sol->instance->Y;
    int *indexPath = sol->indexPath;

    int index1 = genRandom(rndState, 1, nNodes); 
    int index2 = genRandom(rndState, 1, nNodes);
    while (fabsf(index1 - index2) < 2) index2 = genRandom(rndState, 1, nNodes);

    double oldArcsCost = 0;
    double newArcsCost = 0;
    double costEdge1 = 0;   // cost edge (index1-1, index1)
    double costEdge2 = 0;   // cost edge (index1, index1+1)
    double costEdge3 = 0;   // cost edge (index2-1, index2)
    double costEdge4 = 0;   // cost edge (index2, index2+1)

    costEdge1 = computeEdgeCost(X[indexPath[index1-1]], Y[indexPath[index1-1]], X[indexPath[index1]], Y[indexPath[index1]], sol->instance);
    costEdge2 = computeEdgeCost(X[indexPath[index1]], Y[indexPath[index1]], X[indexPath[index1+1]], Y[indexPath[index1+1]], sol->instance);
    costEdge3 = computeEdgeCost(X[indexPath[index2-1]], Y[indexPath[index2-1]], X[indexPath[index2]], Y[indexPath[index2]], sol->instance);
    costEdge4 = computeEdgeCost(X[indexPath[index2]], Y[indexPath[index2]], X[indexPath[index2+1]], Y[indexPath[index2+1]], sol->instance);
    oldArcsCost = costEdge1 + costEdge2 + costEdge3 + costEdge4;

    // now we compute the weights as if we swapped the nodes in the path
    costEdge1 = computeEdgeCost(X[indexPath[index1-1]], Y[indexPath[index1-1]], X[indexPath[index2]], Y[indexPath[index2]], sol->instance);
    costEdge2 = computeEdgeCost(X[indexPath[index2]], Y[indexPath[index2]], X[indexPath[index1+1]], Y[indexPath[index1+1]], sol->instance);
    costEdge3 = computeEdgeCost(X[indexPath[index2-1]], Y[indexPath[index2-1]], X[indexPath[index1]], Y[indexPath[index1]], sol->instance);
    costEdge4 = computeEdgeCost(X[indexPath[index1]], Y[indexPath[index1]], X[indexPath[index2+1]], Y[indexPath[index2+1]], sol->instance);
    newArcsCost = costEdge1 + costEdge2 + costEdge3 + costEdge4;

    SwapInformation swapInfo = {.offset = newArcsCost - oldArcsCost, .index1 = index1, .index2 = index2};

    return swapInfo;
}

static inline void updateTemperature(float *temperature)
{
    *temperature = 0.99 * (*temperature);
}

static inline void normalizeDelta(double *offset, Solution *sol)
{
    *offset = *offset / (sol->cost / sol->instance->nNodes);
}

static inline bool keepMove(double threshold)
{
    if(threshold*RAND_MAX > rand()) return true;
    else return false;
}