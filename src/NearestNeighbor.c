#include "Tsp.h"


#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


typedef struct
{
    Solution *sol;
    
    enum NNFirstNodeOptions startOption;

    pthread_mutex_t getStartNodeMutex;
    pthread_mutex_t saveSolutionMutex;
    size_t startingNode;

    double tlim;
    int nThreads;

}NNThreadsSharedData;

typedef struct
{
    NNThreadsSharedData *thShared;
    unsigned int rndState;
    unsigned long iterCount;
} NNThreadsDataRnd;


static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2);

// finds the closest unvisited node (pathCost is also updated in this method)
static inline size_t findSuccessorID(Solution *sol, size_t iterNum, double *pathCost, unsigned int *state);

static void applyNearestNeighbor(Solution *sol, size_t startingNodeIndex, unsigned int *state);

static void updateBestSolutionNN(Solution *bestSol, Solution *newBest);


static void *nnThread(void *arg);


Solution NearestNeighbor(Instance *inst, enum NNFirstNodeOptions startOption, double timeLimit, bool useThreads)
{
    Solution sol = newSolution(inst);

    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    NNThreadsSharedData th = {
        .sol = &sol,
        .startOption=startOption,
        .startingNode = 0,
        .tlim = timeLimit + (double)start.tv_sec + (double)start.tv_nsec / 1000000000.0,
        .nThreads = inst->params.nThreads
    };

    if (!useThreads)
        th.nThreads = 1;

    unsigned long iterCount = 0;

    if (th.nThreads > 1)
    {
        pthread_mutex_init(&th.getStartNodeMutex, NULL);
        pthread_mutex_init(&th.saveSolutionMutex, NULL);

        NNThreadsDataRnd data[MAX_THREADS];
        pthread_t threads[MAX_THREADS];
        for (int i = 0; i < th.nThreads; i++)
        {
            data[i].thShared = &th; data[i].rndState = rand(); data[i].iterCount = 0;
            pthread_create(&threads[i], NULL, nnThread, (void *)&data);
            LOG(LOG_LVL_DEBUG, "Nearest Neighbour : Thread %d CREATED", i);
        }

        for (int i = 0; i < inst->params.nThreads; i++)
        {
            pthread_join(threads[i], NULL);
            LOG(LOG_LVL_DEBUG, "Nearest Neighbour : Thread %d finished", i);
            iterCount += data[i].iterCount;
        }

        pthread_mutex_destroy(&th.getStartNodeMutex);
        pthread_mutex_destroy(&th.saveSolutionMutex);
    }
    else
    {
        NNThreadsDataRnd data = { .thShared = &th, .rndState = rand(), .iterCount = 0 };
        nnThread(&data);
        iterCount += data.iterCount;
    }

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    sol.execTime = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/sol.execTime);

    return sol;
}

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2)
{
    register float tempf;
    swapElems(sol->X[pos1], sol->X[pos2], tempf);
    swapElems(sol->Y[pos1], sol->Y[pos2], tempf);
    register int tempi;
    swapElems(sol->indexPath[pos1], sol->indexPath[pos2], tempi);
}

static inline size_t findSuccessorID(Solution *sol, size_t iterNum, double *pathCost, unsigned int *state)
{
    Instance *inst = sol->instance;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

    float minVecStore[AVX_VEC_SIZE];
    int minIDsVecStore[AVX_VEC_SIZE];

    // to keep track of the first best distance
    __m256 minVec = _mm256_set1_ps(INFINITY);
    __m256i minIDsVec = _mm256_set1_epi32(-1);

    __m256i ids = _mm256_setr_epi32(0,1,2,3,4,5,6,7), increment = _mm256_set1_epi32(AVX_VEC_SIZE);
    __m256 x1 = _mm256_set1_ps(sol->X[iterNum-1]), y1 = _mm256_set1_ps(sol->Y[iterNum-1]);
    
    // set ids to right starting position
    ids = _mm256_add_epi32(ids, _mm256_set1_epi32((int)iterNum));

    // check if we are working with rounded weights
    for (size_t i = iterNum; i < inst->nNodes; i += AVX_VEC_SIZE)
    {
        // get distance first
        __m256 x2 = _mm256_loadu_ps(&sol->X[i]), y2 = _mm256_loadu_ps(&sol->Y[i]);
        __m256 dist = computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag);

        // get first minimum
        __m256i mask = _mm256_castps_si256(_mm256_cmp_ps(dist, minVec, _CMP_LT_OQ)); // set for "non-signaling" -> if one element is NaN than set as false
        minIDsVec = _mm256_blendv_epi8(minIDsVec, ids, mask);
        minVec = _mm256_min_ps(minVec, dist);

        // increment ids
        ids = _mm256_add_epi32(ids, increment);
    }

    // store vectors
    _mm256_storeu_ps(minVecStore, minVec);
    _mm256_storeu_si256((__m256i*)minIDsVecStore, minIDsVec);

    // sort arrays in ascendant order (size is fixed to 8 so this computation count as linear)
    for (size_t i = 0; i < AVX_VEC_SIZE-1; i++)
    {
        for (size_t j = i+1; j < AVX_VEC_SIZE; j++)
        {
            if (minVecStore[j] < minVecStore[i])
            {
                register float tempf;
                swapElems(minVecStore[i], minVecStore[j], tempf);
                register int tempi;
                swapElems(minIDsVecStore[i], minIDsVecStore[j], tempi);
            }
        }
    }

    // choose successor
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);
    size_t successor = minIDsVecStore[0];
    float pathCostIncrement = minVecStore[0];

    if ((inst->params.graspType != GRASP_NONE) && (rand_r(state) < graspThreshold) && (inst->nNodes - iterNum > 8))
    {
        switch (inst->params.graspType)
        {
        case GRASP_NONE: break; // Prevents warning in compilation

        case GRASP_RANDOM:
            successor = iterNum + (rand_r(state) % (inst->nNodes - iterNum));
            pathCostIncrement = computeEdgeCost(sol->X[iterNum - 1], sol->Y[iterNum - 1], sol->X[successor], sol->Y[successor], ewt, roundFlag);
            break;

        case GRASP_ALMOSTBEST:
        {
            size_t rndIdx = 1;
            for (; rndIdx < AVX_VEC_SIZE-1; rndIdx++)
                if (rand_r(state) > RAND_MAX / 2)
                    break;

            successor = (size_t)minIDsVecStore[rndIdx];
            pathCostIncrement = minVecStore[rndIdx];
            break;
        }
        }
    }

    *pathCost += pathCostIncrement;
    return successor;
}

static void applyNearestNeighbor(Solution *sol, size_t startingNodeIndex, unsigned int *state)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

    // set first node
    swapElementsInSolution(sol, 0, startingNodeIndex);

    // reset cost
    sol->cost = 0.F;

    for (size_t i = 1; i < n-1; i++)
    {
        size_t successorIndex = findSuccessorID(sol, i, &sol->cost, state);

        if (successorIndex >= n)
            throwError(inst, sol, "applyNearestNeighbor: Value of successor isn't applicable: %lu (startNode=%lu)", successorIndex, startingNodeIndex);
        
        // update solution
        swapElementsInSolution(sol, i, successorIndex);
    }

    // add cost of the second to last edge
    sol->cost += computeEdgeCost(sol->X[n-2], sol->Y[n-2], sol->X[n-1], sol->Y[n-1], ewt, roundFlag);

    // add last edge (close the tour)
    sol->X[n] = sol->X[0];
    sol->Y[n] = sol->Y[0];
    sol->indexPath[n] = sol->indexPath[0];
    sol->cost += computeEdgeCost(sol->X[n-1], sol->Y[n-1], sol->X[n], sol->Y[n], ewt, roundFlag);
}

static void updateBestSolutionNN(Solution *bestSol, Solution *newBest)
{
    // check solution when debugging
    if (bestSol->instance->params.logLevel >= LOG_LVL_EVERYTHING)
        checkSolution(newBest);

    LOG(LOG_LVL_LOG, "Found better solution: New cost: %f   Old cost: %f", newBest->cost, bestSol->cost);

    bestSol->cost = newBest->cost;

    { // swap pointers with macro
        register float *temp;
        swapElems(bestSol->X, newBest->X, temp);
        swapElems(bestSol->Y, newBest->Y, temp);
    }
    {
        register int *temp;
        swapElems(bestSol->indexPath, newBest->indexPath, temp);
    }
}




static void *nnThread(void *arg)
{
    NNThreadsDataRnd *thSpecific = (NNThreadsDataRnd*)arg;
    NNThreadsSharedData *thShared = thSpecific->thShared;
    unsigned int state = thSpecific->rndState;
    Solution *sol = thShared->sol;
    Instance *inst = sol->instance;

    // Allocate memory to contain the work-in-progress solution
    Solution tempSol = newSolution(inst);

    double t;
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    t = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;

    // check startingNode even before locking the mutex to avoid useless mutex calls
    while((thShared->startingNode < inst->nNodes) && (t < thShared->tlim))
    {
        size_t iterNode;
        if (thShared->startOption == NN_FIRST_TRYALL)
        {
            pthread_mutex_lock(&thShared->getStartNodeMutex);

            if (thShared->startingNode >= inst->nNodes)
            {
                pthread_mutex_unlock(&thShared->getStartNodeMutex);
                break;
            }

            iterNode = thShared->startingNode;
            thShared->startingNode++;

            pthread_mutex_unlock(&thShared->getStartNodeMutex);
        }
        else
            iterNode = (size_t)rand_r(&state) % inst->nNodes;

        // copy all the nodes in the original order from instance into tempSol.X/Y (Take advantage of the way the memory is allocated for the best and current sol)
        for (size_t i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
            tempSol.X[i] = inst->X[i];

        for (int i = 0; i < (int)inst->nNodes + 1; i++)
            tempSol.indexPath[i] = i;
        
        applyNearestNeighbor(&tempSol, iterNode, &state);

        // check cost first to avoid excessive amount of mutex calls
        if (sol->cost > tempSol.cost)
        {
            pthread_mutex_lock(&thShared->saveSolutionMutex);

            // recheck while having lock on mutex (for syncronization purposes)
            if (sol->cost > tempSol.cost)
                updateBestSolutionNN(thShared->sol, &tempSol);
            
            pthread_mutex_unlock(&thShared->saveSolutionMutex);
        }

        thSpecific->iterCount++;

        //struct timespec currT;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        t = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    }

    destroySolution(&tempSol);

    // avoid closing thread if working in single thread mode
    if (thShared->nThreads == 1)
        return NULL;
    
    pthread_exit(NULL);
}


