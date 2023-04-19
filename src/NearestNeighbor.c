#include "NearestNeighbor.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


typedef struct
{
    Solution *sol;
    pthread_mutex_t nodeLock;
    pthread_mutex_t saveLock;
    size_t startingNode;
}ThreadsData;




static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2);

// finds the closest unvisited node (pathCost is also updated in this method)
static inline size_t findSuccessorID(Solution *sol, size_t iterNum, double *pathCost);

static void applyNearestNeighbor(Solution *sol, size_t startingNodeIndex);

static void updateBestSolutionNN(Solution *bestSol, Solution *newBest);

static void * nnTryAllThread(void *arg);


Solution NearestNeighbor(Instance *inst)
{
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    Solution sol = newSolution(inst);

    // we initialize the seed if it has been passed as argument
    if(inst->params.randomSeed != -1) srand(inst->params.randomSeed);

    // we create and initialize the threaded instance
    ThreadsData th = {.sol = &sol, .startingNode = 0};

    pthread_mutex_init(&th.nodeLock, NULL);
    pthread_mutex_init(&th.saveLock, NULL);
    
    pthread_t threads[inst->params.nThreads];
    for(int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_create(&threads[i], NULL, nnTryAllThread, &th);
        LOG(LOG_LVL_DEBUG, "Nearest Neighbour : Thread %d CREATED", i);
    }

    for(int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        LOG(LOG_LVL_DEBUG, "Nearest Neighbour : Thread %d finished", i);
    }

    pthread_mutex_destroy(&th.nodeLock);
    pthread_mutex_destroy(&th.saveLock);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    sol.execTime = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);

    return sol;
}

Solution NearestNeighborGrasp(Instance *inst, double timeLimit, enum NNThreadsOption option)
{

}

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2)
{
    register float tempf;
    swapElems(sol->X[pos1], sol->X[pos2], tempf);
    swapElems(sol->Y[pos1], sol->Y[pos2], tempf);
    register int tempi;
    swapElems(sol->indexPath[pos1], sol->indexPath[pos2], tempi);
}

static inline size_t findSuccessorID(Solution *sol, size_t iterNum, double *pathCost)
{
    Instance *inst = sol->instance;
    enum edgeWeightType ewt = inst->params.edgeWeightType;
    int roundFlag = inst->params.roundWeightsFlag;

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
                register float tempf = minVecStore[i];
                minVecStore[i] = minVecStore[j];
                minVecStore[j] = tempf;
                register int tempi = minIDsVecStore[i];
                minIDsVecStore[i] = minIDsVecStore[j];
                minIDsVecStore[j] = tempi;
            }
        }
    }

    // choose successor
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);
    size_t successor = minIDsVecStore[0];
    switch (inst->params.graspType)
    {
    case GRASP_NONE:
        *pathCost += minVecStore[0];
        break;

    case GRASP_RANDOM:
        if (rand() < graspThreshold) // wastes an iteration (could be placed outside all as a "big" if however that reduces readability)
        {
            successor = iterNum + (rand() % (inst->nNodes - iterNum));
            *pathCost += computeEdgeCost(sol->X[iterNum-1], sol->Y[iterNum-1], sol->X[successor], sol->Y[successor], ewt, roundFlag);
        }
        else
            *pathCost += minVecStore[0];
        break;
    
    case GRASP_ALMOSTBEST:
        if (rand() < graspThreshold)
        {
            size_t succID = rand() % (AVX_VEC_SIZE - 1);
            if (minIDsVecStore[succID] == -1)
                succID = 0;
                
            successor = minIDsVecStore[succID];
            *pathCost += minVecStore[succID];
        }
        else
            *pathCost += minVecStore[0];
        break;
    }

    return successor;
}

static void applyNearestNeighbor(Solution *sol, size_t startingNodeIndex)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    enum edgeWeightType ewt = inst->params.edgeWeightType;
    int roundFlag = inst->params.roundWeightsFlag;

    // set first node
    swapElementsInSolution(sol, 0, startingNodeIndex);

    // reset cost
    sol->bestCost = 0.F;

    for (size_t i = 1; i < n-1; i++)
    {
        size_t successorIndex = findSuccessorID(sol, i, &sol->bestCost);

        if (successorIndex >= n)
            throwError(inst, sol, "applyNearestNeighbor: Value of successor isn't applicable: %lu (startNode=%lu)", successorIndex, startingNodeIndex);
        
        // update solution
        swapElementsInSolution(sol, i, successorIndex);
    }

    // add cost of the second to last edge
    sol->bestCost += computeEdgeCost(sol->X[n-2], sol->Y[n-2], sol->X[n-1], sol->Y[n-1], ewt, roundFlag);

    // add last edge (close the tour)
    sol->X[n] = sol->X[0];
    sol->Y[n] = sol->Y[0];
    sol->indexPath[n] = sol->indexPath[0];
    sol->bestCost += computeEdgeCost(sol->X[n-1], sol->Y[n-1], sol->X[n], sol->Y[n], ewt, roundFlag);
}

static void updateBestSolutionNN(Solution *bestSol, Solution *newBest)
{
    // check solution when debugging
    if (LOG_LEVEL >= LOG_LVL_EVERYTHING)
        checkSolution(newBest);

    bestSol->bestCost = newBest->bestCost;

    { // swap pointers with macro
        register float *temp;
        swapElems(bestSol->X, newBest->X, temp);
        swapElems(bestSol->Y, newBest->Y, temp);
    }
    {
        register int *temp;
        swapElems(bestSol->indexPath, newBest->indexPath, temp);
    }

    LOG(LOG_LVL_EVERYTHING, "Found better solution starting from node %d, cost: %f", bestSol->indexPath[0], bestSol->bestCost);
}

static void * nnTryAllThread(void *arg)
{
    ThreadsData *th = (ThreadsData *)arg;
    Solution *sol = th->sol;
    Instance *inst = sol->instance;

    // Allocate memory to contain the work-in-progress solution
    Solution tempSol = newSolution(inst);

    // We want the threads to repeat the computation for every node
    // For this we use a mutex on startingNode until it reaches nNodes
    while((pthread_mutex_lock(&th->nodeLock) == 0) && (th->startingNode < inst->nNodes))
    {
        size_t iterNode = th->startingNode;
        th->startingNode++;

        pthread_mutex_unlock(&th->nodeLock);

        // copy all the nodes in the original order from instance into tempSol.X/Y (Take advantage of the way the memory is allocated for the best and current sol)
        for (size_t i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
            tempSol.X[i] = inst->X[i];

        for (int i = 0; i < (int)inst->nNodes + 1; i++)
            tempSol.indexPath[i] = i;
        
        applyNearestNeighbor(&tempSol, iterNode);

        // check bestCost first to avoid excessive amount of mutex calls
        if (sol->bestCost > tempSol.bestCost)
        {
            pthread_mutex_lock(&th->saveLock);

            // recheck while having lock on mutex (for syncronization purposes)
            if (sol->bestCost > tempSol.bestCost)
                updateBestSolutionNN(th->sol, &tempSol);
            
            pthread_mutex_unlock(&th->saveLock);
        }
    }
    
    // mutex is locked in while condition
    pthread_mutex_unlock(&th->nodeLock);

    destroySolution(&tempSol);

    pthread_exit(NULL);
}



