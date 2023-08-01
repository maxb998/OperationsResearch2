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

    if (startOption == NN_FIRST_TRYALL) // syncronization in the selection of the starting node of the solution is only necessary when tryall option is used
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

    if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
    {   // this memory allocation is useful, but not strictly necessary when using avx (if we don't want to use it switch to the gather instructions)
        thSpecific.X = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float));
        if (!thSpecific.X)
            throwError(inst, &thSpecific.workingSol, "NearestNeighbor -> initThreadSpecificData: Failed to allocate memory");
        thSpecific.Y = &thSpecific.X[inst->nNodes + AVX_VEC_SIZE];
    }
    else
        thSpecific.X = thSpecific.Y = NULL;

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    free(thSpecific->X);
    thSpecific->X = thSpecific->Y = NULL;
}

typedef struct 
{
    int node;
    float cost;
} SuccessorData;


// Executes Nearest Neighbor algorithm multiple times until time limit (designed to work in parallel)
static void *loopNearestNeighbor(void *arg);

Solution NearestNeighbor(Instance *inst, enum NNFirstNodeOptions startOption, double timeLimit, int nThreads)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // check input
    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError(inst, NULL, "NearestNeighbor: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = inst->params.nThreads;

    // Create data structures and start threads
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

        // check cost before locking mutex to avoid excessive amount of mutex calls
        if (thShared->bestSol.cost > thSpecific->workingSol.cost)
        {
            pthread_mutex_lock(&thShared->saveSolutionMutex);

            // recheck while having lock on mutex (for syncronization purposes)
            if (thShared->bestSol.cost > thSpecific->workingSol.cost)
                updateBestSolution(thSpecific);
            
            pthread_mutex_unlock(&thShared->saveSolutionMutex);
        }

        thSpecific->iterCount++;

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

// finds the closest unvisited node using vectorized(SIMD) instructions
static inline SuccessorData findSuccessorVectorized(ThreadSpecificData *thSpecific, int lastPos);

// finds the closest unvisited node using normal(SISD) instructions
static inline SuccessorData findSuccessorBase(ThreadSpecificData *thSpecific, int lastAddedPos);

static void applyNearestNeighbor(ThreadSpecificData *thSpecific, int firstNode)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    int *workSolPath = thSpecific->workingSol.indexPath;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

    // Reset thSpecific->[X,Y,workingSol]
    if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
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
        SuccessorData successor;

        // if GRASP_RANDOM is used and this iteration should return a successor selected completely at random then the we just skip the computation imposed by "findSuccessor" functions
        if ((inst->params.graspType == GRASP_RANDOM) && (rand_r(&thSpecific->rndState) < graspThreshold))
        {
            successor.node = genRandom(&thSpecific->rndState, (i+1), n);
            switch (COMPUTATION_TYPE)
            {
            case COMPUTE_OPTION_AVX:
                successor.cost = computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[successor.node], thSpecific->Y[successor.node], ewt, roundFlag);
                break;
            
            case COMPUTE_OPTION_BASE:
                successor.cost = computeEdgeCost(inst->X[workSolPath[i]], inst->Y[workSolPath[i]], inst->X[workSolPath[successor.node]], inst->Y[workSolPath[successor.node]], ewt, roundFlag);
                break;

            case COMPUTE_OPTION_USE_COST_MATRIX:
                successor.cost = inst->edgeCostMat[workSolPath[i] * n + successor.node];
                break;
            }
        }
        else
        {
            switch (COMPUTATION_TYPE)
            {
            case COMPUTE_OPTION_AVX:
                successor = findSuccessorVectorized(thSpecific, i);
                break;
            
            case COMPUTE_OPTION_BASE:
            case COMPUTE_OPTION_USE_COST_MATRIX:
                successor = findSuccessorBase(thSpecific, i);
                break;
            }
        }

        // simple debugging check. can be removed, but saved a lot of headaches so it's going to stay there
        if ((successor.node >= n) || (successor.node < 0))
            throwError(inst, &thSpecific->workingSol, "applyNearestNeighbor: Value of successor isn't applicable: %d (startNode=%d)", successor.node, firstNode);

        // update solution
        swapElemsInThSpecific(thSpecific, i+1, successor.node);
        thSpecific->workingSol.cost += successor.cost;
    }

    // add cost of the two remaining edges
    switch (COMPUTATION_TYPE)
    {
    case COMPUTE_OPTION_AVX:
        thSpecific->workingSol.cost += computeEdgeCost(thSpecific->X[n-2], thSpecific->Y[n-2], thSpecific->X[n-1], thSpecific->Y[n-1], ewt, roundFlag) +
                                       computeEdgeCost(thSpecific->X[n-1], thSpecific->Y[n-1], thSpecific->X[0],   thSpecific->Y[0],   ewt, roundFlag);
        break;

    case COMPUTE_OPTION_BASE:
        thSpecific->workingSol.cost += computeEdgeCost(inst->X[workSolPath[n-2]], inst->Y[workSolPath[n-2]], inst->X[workSolPath[n-1]], inst->Y[workSolPath[n-1]], ewt, roundFlag) + 
                                       computeEdgeCost(inst->X[workSolPath[n-1]], inst->Y[workSolPath[n-1]], inst->X[workSolPath[0]],   inst->Y[workSolPath[0]],   ewt, roundFlag);
        break;

    case COMPUTE_OPTION_USE_COST_MATRIX:
        thSpecific->workingSol.cost += inst->edgeCostMat[workSolPath[n-2] * n + workSolPath[n-1]] + 
                                       inst->edgeCostMat[workSolPath[n-1] * n + workSolPath[0]];
        break;
    }
    
}

static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2)
{
    if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX) //thSpecific->X in not used otherwise
    {
        register float tempFloat;
        swapElems(thSpecific->X[pos1], thSpecific->X[pos2], tempFloat);
        swapElems(thSpecific->Y[pos1], thSpecific->Y[pos2], tempFloat);
    }
    register int tempInt;
    swapElems(thSpecific->workingSol.indexPath[pos1], thSpecific->workingSol.indexPath[pos2], tempInt);
}

static void updateBestSolution(ThreadSpecificData *thSpecific)
{
    Solution *bestSol = &thSpecific->thShared->bestSol;
    Solution *newBest = &thSpecific->workingSol;

    // update cost as first things
    bestSol->cost = newBest->cost;

    // check solution when debugging
    if (bestSol->instance->params.logLevel >= LOG_LVL_DEBUG)
        if (!checkSolution(newBest))
		    throwError(newBest->instance, newBest, "updateBestSolution: newBest Solution is not valid");

    LOG(LOG_LVL_LOG, "Found better solution: cost = %f", newBest->cost);

    register int *temp;
    swapElems(bestSol->indexPath, newBest->indexPath, temp);
}

static inline SuccessorData findSuccessorVectorized(ThreadSpecificData *thSpecific, int lastAddedPos)
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
        __m256 dist = computeSquaredEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag);

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

    SuccessorData succ = { .cost=sqrtf(minVecStore[successorSubIndex]), .node=minIDsVecStore[successorSubIndex] };

    return succ;
}

static inline SuccessorData findSuccessorBase(ThreadSpecificData *thSpecific, int lastAddedPos)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    int *indexPath = thSpecific->workingSol.indexPath;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

    // initialize array with stored best successors(1st best, 2nd best, 3rd best, ...)
    SuccessorData bestSuccs[BASE_GRASP_BEST_SAVE_BUFFER_SIZE];
    for (int i = 0; i < BASE_GRASP_BEST_SAVE_BUFFER_SIZE; i++)
    {
        bestSuccs[i].cost = INFINITY;
        bestSuccs[i].node = -1;
    }
    
    // hope the compiler loads into the registers these variables that are accessed every time in the loop
    int lastAddedIndex = indexPath[lastAddedPos];
    float x1, y1;
    if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
    {
        x1 = inst->X[indexPath[lastAddedPos]];
        y1 = inst->Y[indexPath[lastAddedPos]];
    }

    int posToAdd = lastAddedPos + 1;

    for (int node = posToAdd ; node < n; node++)
    {   
        SuccessorData currentSucc = { .node=node };

        if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            currentSucc.cost = computeSquaredEdgeCost(x1, y1, inst->X[indexPath[currentSucc.node]], inst->Y[indexPath[currentSucc.node]], ewt, roundFlag);
        else
            currentSucc.cost = inst->edgeCostMat[lastAddedIndex * n + indexPath[currentSucc.node]];

        for (int i = 0; i < BASE_GRASP_BEST_SAVE_BUFFER_SIZE; i++)
        {
            if (currentSucc.cost < bestSuccs[i].cost)
            {
                SuccessorData temp;
                swapElems(bestSuccs[i], currentSucc, temp);
            }
        }
    }

    // choose successor
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);
    int successorSubIndex = 0;

    // if using GRASP_ALMOSTBEST this is the time to select with some probability one of the nodes saved in bestSuccs[]
    if ((inst->params.graspType == GRASP_ALMOSTBEST) && (rand_r(&thSpecific->rndState) < graspThreshold) && (n - posToAdd > BASE_GRASP_BEST_SAVE_BUFFER_SIZE))
    {
        successorSubIndex = 1;
        for (; successorSubIndex < BASE_GRASP_BEST_SAVE_BUFFER_SIZE - 1; successorSubIndex++)
            if (rand_r(&thSpecific->rndState) > RAND_MAX / 2)
                break;
    }

    if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        bestSuccs[successorSubIndex].cost = sqrtf(bestSuccs[successorSubIndex].cost);

    return bestSuccs[successorSubIndex];
}
