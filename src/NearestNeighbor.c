#include "Tsp.h"


#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


typedef struct
{
    Solution bestSol;
    
    enum NNFirstNodeOptions startOption;

    pthread_mutex_t getStartNodeMutex;
    pthread_mutex_t saveSolutionMutex;
    int startingNode;

    double timeLimit;

} ThreadSharedData;

static ThreadSharedData initThreadSharedData (Instance *inst, enum NNFirstNodeOptions startOption, double timeLimit)
{
    ThreadSharedData thShared = {
        .startOption = startOption,
        .startingNode = 0,
        .timeLimit = timeLimit
    };

    if (startOption == NN_FIRST_TRYALL)
        if (pthread_mutex_init(&thShared.getStartNodeMutex, NULL)) throwError(inst, NULL, "NearestNeighbor -> initThreadSharedData: Failed to initialize getStartNodeMutex");
    if (pthread_mutex_init(&thShared.saveSolutionMutex, NULL)) throwError(inst, NULL, "NearestNeighbor -> initThreadSharedData: Failed to initialize saveSolutionMutex");

    thShared.bestSol = newSolution(inst);

    return thShared;
}

static void destroyThreadSharedData (ThreadSharedData *thShared)
{
    if (thShared->startOption == NN_FIRST_TRYALL)
        if (pthread_mutex_init(&thShared->getStartNodeMutex, NULL)) throwError(thShared->bestSol.instance, &thShared->bestSol, "NearestNeighbor -> destroyThreadSharedData: Failed to destroy getStartNodeMutex");
    if (pthread_mutex_init(&thShared->saveSolutionMutex, NULL)) throwError(thShared->bestSol.instance, &thShared->bestSol, "NearestNeighbor -> destroyThreadSharedData: Failed to destroy saveSolutionMutex");

    destroySolution(&thShared->bestSol);
}

typedef struct
{
    ThreadSharedData *thShared;
    unsigned int rndState;
    int iterCount;

    float *X;
    float *Y;
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

    thSpecific.X = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float));
    if (!thSpecific.X)
        throwError(inst, &thSpecific.workingSol, "NearestNeighbor -> initThreadSpecificData: Failed to allocate memory");
    thSpecific.Y = &thSpecific.X[inst->nNodes + AVX_VEC_SIZE];

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    free(thSpecific->X);
    thSpecific->X = thSpecific->Y = NULL;
}

// Executes Nearest Neighbor algorithm multiple times until time limit (designed to work in parallel)
static void *loopNearestNeighbor(void *arg);

Solution NearestNeighbor(Instance *inst, enum NNFirstNodeOptions startOption, double timeLimit, int nThreads)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError(inst, NULL, "NearestNeighbor: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = inst->params.nThreads;

    ThreadSharedData thShared = initThreadSharedData(inst, startOption, startTime + timeLimit);

    ThreadSpecificData thSpecifics[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, rand());
        pthread_create(&threads[i], NULL, loopNearestNeighbor, (void *)&thSpecifics[i]);
    }

    int iterCount = 0;
    for (int i = 0; i < nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        iterCount += thSpecifics[i].iterCount;
        destroyThreadSpecificData(&thSpecifics[i]);
    }

    Solution outputSol = newSolution(inst);
    cloneSolution(&thShared.bestSol, &outputSol);

    destroyThreadSharedData(&thShared);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    outputSol.execTime = cvtTimespec2Double(timeStruct) - startTime;

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/outputSol.execTime);

    return outputSol;
}

// select node to use as starting point for a Nearest Neighbor iteration depending on startOption
static inline int getStartingNode(ThreadSpecificData *thSpecific);

// Apply Nearest Neighbor algorithm to a solution that has the starting node in the first position and all others next
static void applyNearestNeighbor(ThreadSpecificData *thSpecific, int firstNode);

// Update thShared->bestSol when a better solution is found
static void updateBestSolution(ThreadSpecificData *thSpecific);

static void *loopNearestNeighbor(void *arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    // check startingNode even before locking the mutex to avoid useless mutex calls
    while((thShared->startingNode < n) && (currentTime < thShared->timeLimit))
    {
        int startNode = getStartingNode(thSpecific);
        
        applyNearestNeighbor(thSpecific, startNode);

        // check cost first to avoid excessive amount of mutex calls
        if (thShared->bestSol.cost > thSpecific->workingSol.cost)
        {
            pthread_mutex_lock(&thShared->saveSolutionMutex);

            // recheck while having lock on mutex (for syncronization purposes)
            if (thShared->bestSol.cost > thSpecific->workingSol.cost)
                updateBestSolution(thSpecific);
            
            pthread_mutex_unlock(&thShared->saveSolutionMutex);
        }

        thSpecific->iterCount++;

        //struct timespec timeStruct;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }
    
    return NULL;
}

static inline int getStartingNode(ThreadSpecificData *thSpecific)
{
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thShared->bestSol.instance;
    int n = inst->nNodes;

    int iterNode;
    if (thShared->startOption == NN_FIRST_TRYALL)
    {
        pthread_mutex_lock(&thShared->getStartNodeMutex);

        if (thShared->startingNode >= n)
        {
            if (inst->params.graspType == GRASP_NONE) // tested all options -> quit
                return -1;
            else // if grasp is being used then start again from first node until time limit is hit
                thShared->startingNode = 0;
        }

        iterNode = thShared->startingNode;
        thShared->startingNode++;

        pthread_mutex_unlock(&thShared->getStartNodeMutex);
    }
    else
        iterNode = genRandom(&thSpecific->rndState, 0, n);
    
    return iterNode;
}

// Swap elements in thSpecific.X, thSpecific.Y and thSpecific.workingSol.indexPath
static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2);

// finds the closest unvisited node (pathCost is also updated in this method)
static inline int findSuccessorID(ThreadSpecificData *thSpecific, int lastPos);

static void applyNearestNeighbor(ThreadSpecificData *thSpecific, int firstNode)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

    // Reset thSpecific->[X,Y,workingSol]
    for (int i = 0; i < (n + AVX_VEC_SIZE) * 2; i++)
        thSpecific->X[i] = inst->X[i];
    for (int i = 0; i < n + AVX_VEC_SIZE; i++)
        thSpecific->workingSol.indexPath[i] = i;

    // set first node
    swapElemsInThSpecific(thSpecific, 0, firstNode);

    // reset cost
    thSpecific->workingSol.cost = 0;

    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);
    for (int i = 0; i < n-2; i++)
    {
        int successor;

        if ((inst->params.graspType == GRASP_RANDOM) && (rand_r(&thSpecific->rndState) < graspThreshold))
        {
            successor = genRandom(&thSpecific->rndState, (i+1), n);
            thSpecific->workingSol.cost += computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[successor], thSpecific->Y[successor], ewt, roundFlag);
        }
        else
            successor = findSuccessorID(thSpecific, i);

        if ((successor >= n) || (successor < 0))
            throwError(inst, &thSpecific->workingSol, "applyNearestNeighbor: Value of successor isn't applicable: %d (startNode=%d)", successor, firstNode);
        
        // update solution
        swapElemsInThSpecific(thSpecific, i+1, successor);
    }

    // add cost of the second to last edge
    thSpecific->workingSol.cost += computeEdgeCost(thSpecific->X[n-2], thSpecific->Y[n-2], thSpecific->X[n-1], thSpecific->Y[n-1], ewt, roundFlag);

    // add last edge (close the tour)
    //thSpecific->workingSol.indexPath[n] = thSpecific->workingSol.indexPath[0];
    thSpecific->workingSol.cost += computeEdgeCost(thSpecific->X[n-1], thSpecific->Y[n-1], thSpecific->X[0], thSpecific->Y[0], ewt, roundFlag);
}

static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2)
{
    register float tempFloat;
    swapElems(thSpecific->X[pos1], thSpecific->X[pos2], tempFloat);
    swapElems(thSpecific->Y[pos1], thSpecific->Y[pos2], tempFloat);
    register int tempInt;
    swapElems(thSpecific->workingSol.indexPath[pos1], thSpecific->workingSol.indexPath[pos2], tempInt);
}

static inline int findSuccessorID(ThreadSpecificData *thSpecific, int lastAddedPos)
{
    Instance *inst = thSpecific->workingSol.instance;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

    // to keep track of the first best distance
    __m256 minVec = _mm256_set1_ps(INFINITY);
    __m256i minIDsVec = _mm256_set1_epi32(-1);

    __m256i ids = _mm256_setr_epi32(0,1,2,3,4,5,6,7), increment = _mm256_set1_epi32(AVX_VEC_SIZE);
    __m256 x1 = _mm256_broadcast_ss(&thSpecific->X[lastAddedPos]), y1 = _mm256_broadcast_ss(&thSpecific->Y[lastAddedPos]);
    
    int posToAdd = lastAddedPos + 1;

    // set ids to right starting position
    ids = _mm256_add_epi32(ids, _mm256_set1_epi32(posToAdd));

    // check if we are working with rounded weights
    for (int i = posToAdd; i < inst->nNodes; i += AVX_VEC_SIZE)
    {
        // get distance first
        __m256 x2 = _mm256_loadu_ps(&thSpecific->X[i]), y2 = _mm256_loadu_ps(&thSpecific->Y[i]);
        __m256 dist = computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag);

        // get first minimum
        __m256i mask = _mm256_castps_si256(_mm256_cmp_ps(dist, minVec, _CMP_LT_OQ)); // set for "non-signaling" -> if one element is NaN than set as false
        minIDsVec = _mm256_blendv_epi8(minIDsVec, ids, mask);
        minVec = _mm256_min_ps(minVec, dist);

        // increment ids
        ids = _mm256_add_epi32(ids, increment);
    }

    // store vectors
    float minVecStore[AVX_VEC_SIZE];
    int minIDsVecStore[AVX_VEC_SIZE];
    _mm256_storeu_ps(minVecStore, minVec);
    _mm256_storeu_si256((__m256i*)minIDsVecStore, minIDsVec);

    // argsort for grasp and get minium(since they are only 8 elements it's fast)
    int sortedArgs[AVX_VEC_SIZE];
    argsort(minVecStore, sortedArgs, AVX_VEC_SIZE);

    // choose successor
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);
    int successorSubIndex = sortedArgs[0];

    if ((inst->params.graspType == GRASP_ALMOSTBEST) && (rand_r(&thSpecific->rndState) < graspThreshold) && (inst->nNodes - posToAdd > AVX_VEC_SIZE))
    {
        int rndIdx = 1;
        for (; rndIdx < AVX_VEC_SIZE - 1; rndIdx++)
            if (rand_r(&thSpecific->rndState) > RAND_MAX / 2)
                break;

        successorSubIndex = sortedArgs[rndIdx];
    }

    thSpecific->workingSol.cost += minVecStore[successorSubIndex];
    return minIDsVecStore[successorSubIndex];
}

static void updateBestSolution(ThreadSpecificData *thSpecific)
{
    Solution *bestSol = &thSpecific->thShared->bestSol;
    Solution *newBest = &thSpecific->workingSol;

    // check solution when debugging
    if (bestSol->instance->params.logLevel >= LOG_LVL_DEBUG)
        if (!checkSolution(newBest))
		    throwError(newBest->instance, newBest, "updateBestSolutionNN: newBest Solution is not valid");

    LOG(LOG_LVL_LOG, "Found better solution: cost = %f", newBest->cost);

    bestSol->cost = newBest->cost;

    register int *temp;
    swapElems(bestSol->indexPath, newBest->indexPath, temp);
}
