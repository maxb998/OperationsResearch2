#include "Tsp.h"


#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


//#define DEBUG

typedef struct
{
    Solution bestSol;
    float bestCost;
    double timeLimit;

    pthread_mutex_t saveSolutionMutex;
    int startingNode;

} ThreadSharedData;

typedef struct
{
    ThreadSharedData *thShared;
    unsigned int rndState;
    int iterCount;

    __uint128_t localBestCost;
    int *localBestPath;
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        float *localBestX;
        float *localBestY;
    #endif

    __uint128_t cost;
    int *path;
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        float *X;
        float *Y;
    #endif

} ThreadSpecificData;

typedef struct 
{
    int node;
    float cost;
} SuccessorData;


// Setup internal variables and initializes mutexes
static ThreadSharedData initThreadSharedData (Instance *inst, double timeLimit);

// Destroy mutexes
static void destroyThreadSharedData (ThreadSharedData *thShared);

// Setup internal variables and allocate memory in thSpecific.X/Y if using COMPUTE_OPTION_AVX
static ThreadSpecificData *initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState);

// Deallocate memory pointed by thSpecific.X/Y if necessary
static void destroyThreadSpecificData(ThreadSpecificData *thSpecific);

// Executes Nearest Neighbor algorithm multiple times until time limit (designed to work in parallel)
static void *loopNearestNeighbor(void *arg);

// select node to use as starting point for a Nearest Neighbor iteration depending on startOption
static inline int getStartingNode(ThreadSpecificData *thSpecific);

// Apply Nearest Neighbor algorithm to a solution that has the starting node in the first position and all others next
static void applyNearestNeighbor(ThreadSpecificData *thSpecific, int firstNode);

// Swap elements in thSpecific.X, thSpecific.Y and thSpecific.workingSol.indexPath
static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2);

#if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
// finds the closest unvisited node using vectorized(SIMD) instructions
static inline SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int lastPos);
#elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
// finds the closest unvisited node using normal(SISD) instructions
static inline SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int lastAddedPos);
#endif


Solution NearestNeighbor(Instance *inst, double timeLimit)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    if (inst->params.graspChance == -1)
    {
        if (inst->params.graspType == GRASP_ALMOSTBEST)
        {
            // y = a x^b.  a and b obtained by using linear regression on a dataset composed of 3264 runs
            inst->params.graspChance = exp(0.68515) * pow(inst->nNodes, -0.6464);
        }
        else
        {
            // y = x^b.   b obtained by using linear regression on a dataset composed of 3264 runs
            inst->params.graspChance = pow(inst->nNodes, -0.87173);
        }
        // cap at 0.5
        if (inst->params.graspChance > 0.5)
            inst->params.graspChance = 0.5;
        LOG(LOG_LVL_NOTICE, "Selected Grasp chance: %lf", inst->params.graspChance);
    }

    // Create data structures and start threads
    ThreadSharedData thShared = initThreadSharedData(inst, startTime + timeLimit);
    ThreadSpecificData **thSpecifics = malloc(inst->params.nThreads * sizeof(ThreadSpecificData*));
    if (thSpecifics == NULL)
        throwError("Failed to allocate memory for **thSpecifics");
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, rand());
        pthread_create(&threads[i], NULL, loopNearestNeighbor, (void *)thSpecifics[i]);
    }

    int iterCount = 0;
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        iterCount += thSpecifics[i]->iterCount;
    }

    Solution outputSol = newSolution(inst);
    cloneSolution(&thShared.bestSol, &outputSol);

    for (int i = 0; i < inst->params.nThreads; i++)
        destroyThreadSpecificData(thSpecifics[i]);
    free(thSpecifics);

    destroyThreadSharedData(&thShared);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    outputSol.execTime = cvtTimespec2Double(timeStruct) - startTime;

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/outputSol.execTime);

    return outputSol;
}


static ThreadSharedData initThreadSharedData (Instance *inst, double timeLimit)
{
    ThreadSharedData thShared = {
        .startingNode = 0,
        .bestCost=INFINITY,
        .timeLimit=timeLimit
    };

    if (pthread_mutex_init(&thShared.saveSolutionMutex, NULL)) throwError("NearestNeighbor -> initThreadSharedData: Failed to initialize saveSolutionMutex");

    thShared.bestSol = newSolution(inst);

    return thShared;
}

static void destroyThreadSharedData (ThreadSharedData *thShared)
{
    if (pthread_mutex_init(&thShared->saveSolutionMutex, NULL)) throwError("NearestNeighbor -> destroyThreadSharedData: Failed to destroy saveSolutionMutex");
}

static ThreadSpecificData *initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState)
{
    Instance *inst = thShared->bestSol.instance;
    size_t memToAlloc = sizeof(ThreadSpecificData) + (inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(int);
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        memToAlloc += (inst->nNodes + AVX_VEC_SIZE) * 4 * sizeof(float);
    #endif

    ThreadSpecificData *thSpecific = malloc(memToAlloc);
    if (thSpecific == NULL)
        throwError("initThreadSpecificData: Failed to allocate memory");

    thSpecific->thShared = thShared;
    thSpecific->rndState = rndState;
    thSpecific->iterCount = 0;
    thSpecific->localBestCost = -1;

    thSpecific->localBestPath = (int*)&thSpecific[1];
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        thSpecific->localBestX = (float*)&thSpecific->localBestPath[inst->nNodes + AVX_VEC_SIZE];
        thSpecific->localBestY = &thSpecific->localBestX[inst->nNodes + AVX_VEC_SIZE];
        thSpecific->path = (int*)&thSpecific->localBestY[inst->nNodes + AVX_VEC_SIZE];
        thSpecific->X = (float*)&thSpecific->path[inst->nNodes + AVX_VEC_SIZE];
        thSpecific->Y = &thSpecific->X[inst->nNodes + AVX_VEC_SIZE];
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        thSpecific->path = &thSpecific->localBestPath[inst->nNodes + AVX_VEC_SIZE];
    #endif

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    free(thSpecific);
}

static void *loopNearestNeighbor(void *arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thShared->bestSol.instance;
    int n = inst->nNodes;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    // set thSpecific->[path,X,Y] and localBest
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        for (int i = 0; i < (n + AVX_VEC_SIZE) * 2; i++)
            thSpecific->X[i] = inst->X[i];
        for (int i = 0; i < (n + AVX_VEC_SIZE) * 2; i++)
            thSpecific->localBestX[i] = inst->X[i];
    #endif
    for (int i = 0; i < n + AVX_VEC_SIZE; i++)
        thSpecific->path[i] = i;
    for (int i = 0; i < n + AVX_VEC_SIZE; i++)
        thSpecific->localBestPath[i] = i;

    // check startingNode even before locking the mutex to avoid useless mutex calls
    while(currentTime < thShared->timeLimit)
    {
        int startNode = getStartingNode(thSpecific);
        if (startNode == -1)
            break;
        
        if (inst->params.nnFirstNodeOption == NN_FIRST_TRYALL)
        {
            // reset needed in case we want a fully deterministic algorithm(starting order inside thSpecific.[X,Y,path] influences the output)
            #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                for (int i = 0; i < (n + AVX_VEC_SIZE) * 2; i++)
                    thSpecific->X[i] = inst->X[i];
            #endif
            for (int i = 0; i < n + AVX_VEC_SIZE; i++)
                thSpecific->path[i] = i;

            // commented code below allows a non-deterministic approach (randomness introduced by not resetting thSpecific.[X,Y,path] to starting position)
            /*int startNodeIndex = 0;
            for (; startNodeIndex < n; startNodeIndex++)
                if (thSpecific->path[startNodeIndex] == startNode)
                    break;
            startNode = startNodeIndex;*/
        }

        applyNearestNeighbor(thSpecific, startNode);

        // check cost before locking mutex to avoid excessive amount of mutex calls
        if (thSpecific->cost < thSpecific->localBestCost)
        {
            thSpecific->localBestCost = thSpecific->cost;
            swapElems(thSpecific->path, thSpecific->localBestPath)
            #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                swapElems(thSpecific->X, thSpecific->localBestX)
                swapElems(thSpecific->Y, thSpecific->localBestY)
            #endif

            if ((float)cvtCost2Double(thSpecific->cost) < thShared->bestCost)
            {
                pthread_mutex_lock(&thShared->saveSolutionMutex);
                if ((float)cvtCost2Double(thSpecific->cost) < thShared->bestCost)
                {
                    thShared->bestCost = (float)cvtCost2Double(thSpecific->cost);
                    thShared->bestSol.indexPath = thSpecific->localBestPath;
                    thShared->bestSol.cost = thSpecific->cost;
                    #ifdef DEBUG
                        if (!checkSolution(&thShared->bestSol))
                            throwError("New solution not correct");
                    #endif
                    LOG(LOG_LVL_INFO, "Found better solution starting from node %d\t with cost: %f", thSpecific->localBestPath[0], thShared->bestCost);
                }
                pthread_mutex_unlock(&thShared->saveSolutionMutex);
            }
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
    if (inst->params.nnFirstNodeOption == NN_FIRST_TRYALL)
    {
        register int startNodeVal = thShared->startingNode;
        if (startNodeVal >= n - 1)
        {
            if (inst->params.graspType == GRASP_NONE) // tested all options -> quit
                return -1;
            else // if grasp is being used then start again from first node until time limit is hit
                startNodeVal = 0;
        }

        iterNode = startNodeVal;
        thShared->startingNode = startNodeVal + 1;
    }
    else
        iterNode = genRandom(&thSpecific->rndState, 1, n);
    
    return iterNode;
}

static void applyNearestNeighbor(ThreadSpecificData *thSpecific, int firstNode)
{
    Instance *inst = thSpecific->thShared->bestSol.instance;
    int n = inst->nNodes;

    // set first node
    swapElemsInThSpecific(thSpecific, 0, firstNode);

    // reset cost
    thSpecific->cost = 0;

    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);
    for (int i = 0; i < n-2; i++)
    {
        SuccessorData successor;

        // if GRASP_RANDOM is used and this iteration should return a successor selected completely at random then the we just skip the computation imposed by "findSuccessor" functions
        if ((inst->params.graspType == GRASP_RANDOM) && (rand_r(&thSpecific->rndState) < graspThreshold))
        {
            successor.node = genRandom(&thSpecific->rndState, (i+1), n);

            #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                successor.cost = computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[successor.node], thSpecific->Y[successor.node], inst);
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                successor.cost = inst->edgeCostMat[(size_t)thSpecific->path[i] * (size_t)n + (size_t)thSpecific->path[successor.node]];
            #endif
        }
        else
        {
            #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
                successor = findSuccessor(thSpecific, i);
            #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
                successor = findSuccessor(thSpecific, i);
            #endif
        }

        // simple debugging check. can be removed, but saved a lot of headaches so it's going to stay there
        #ifdef DEBUG
            if ((successor.node >= n) || (successor.node < 0))
                throwError("applyNearestNeighbor: Value of successor isn't applicable: %d (startNode=%d)", successor.node, firstNode);
        #endif

        // update solution
        swapElemsInThSpecific(thSpecific, i+1, successor.node);
        thSpecific->cost += cvtFloat2Cost(successor.cost);
    }

    float secondToLastCost, lastCost;

    // add cost of the two remaining edges
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        secondToLastCost = computeEdgeCost(thSpecific->X[n - 2], thSpecific->Y[n - 2], thSpecific->X[n - 1], thSpecific->Y[n - 1], inst);
        lastCost = computeEdgeCost(thSpecific->X[n - 1], thSpecific->Y[n - 1], thSpecific->X[0], thSpecific->Y[0], inst);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        secondToLastCost = inst->edgeCostMat[(size_t)thSpecific->path[n-2] * (size_t)n + (size_t)thSpecific->path[n-1]];
        lastCost = inst->edgeCostMat[(size_t)thSpecific->path[n-1] * (size_t)n + (size_t)thSpecific->path[0]];
    #endif

    thSpecific->cost += cvtFloat2Cost(secondToLastCost);
    thSpecific->cost += cvtFloat2Cost(lastCost);
}

static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2)
{
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)) //thSpecific->X in not used otherwise
        swapElems(thSpecific->X[pos1], thSpecific->X[pos2])
        swapElems(thSpecific->Y[pos1], thSpecific->Y[pos2])
    #endif
    swapElems(thSpecific->path[pos1], thSpecific->path[pos2])
}

#if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
static inline SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int lastAddedPos)
{
    Instance *inst = thSpecific->thShared->bestSol.instance;

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
        __m256 dist = noSquaredRootEdgeCost_VEC(x1, y1, x2, y2, inst);

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

    // choose successor
    int minIndex = 0;
    for (int i = 1; i < AVX_VEC_SIZE; i++)
        if (minVecStore[i] < minVecStore[minIndex])
            minIndex = i;

    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);
    if ((inst->params.graspType == GRASP_ALMOSTBEST) && (rand_r(&thSpecific->rndState) < graspThreshold) && (inst->nNodes - posToAdd > AVX_VEC_SIZE))
    {
        minVecStore[minIndex] = INFINITY;
        for (int counter = 1; counter < AVX_VEC_SIZE; counter++)
        {
            if (rand_r(&thSpecific->rndState) < graspThreshold)
            {
                for (int i = 0; i < AVX_VEC_SIZE; i++)
                    if (minVecStore[i] < minVecStore[minIndex])
                        minIndex = i;
                minVecStore[minIndex] = INFINITY;
            }
        }
    }

    int nextPos = minIDsVecStore[minIndex];
    SuccessorData succ = {
        .cost=computeEdgeCost(thSpecific->X[lastAddedPos], thSpecific->Y[lastAddedPos], thSpecific->X[nextPos], thSpecific->Y[nextPos], inst), 
        .node=nextPos
    };

    return succ;
}
#elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
static inline SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int lastAddedPos)
{
    Instance *inst = thSpecific->thShared->bestSol.instance;
    int n = inst->nNodes;

    // initialize array with stored best successors(1st best, 2nd best, 3rd best, ...)
    SuccessorData bestSuccs[BASE_GRASP_BEST_SAVE_BUFFER_SIZE];
    for (int i = 0; i < BASE_GRASP_BEST_SAVE_BUFFER_SIZE; i++)
    {
        bestSuccs[i].cost = INFINITY;
        bestSuccs[i].node = -1;
    }
    
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
    int lastAddedIndex = thSpecific->path[lastAddedPos];
    #endif

    int posToAdd = lastAddedPos + 1;

    for (int node = posToAdd; node < n; node++)
    {   
        SuccessorData currentSucc = { .node=node };

        #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            currentSucc.cost = noSquaredRootEdgeCost(thSpecific->X[lastAddedPos], thSpecific->Y[lastAddedPos], thSpecific->X[node], thSpecific->Y[node], inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            currentSucc.cost = inst->edgeCostMat[(size_t)lastAddedIndex * (size_t)n + (size_t)thSpecific->path[node]];
        #endif

        for (int i = 0; i < BASE_GRASP_BEST_SAVE_BUFFER_SIZE; i++)
        {
            if (currentSucc.cost < bestSuccs[i].cost)
                swapElems(bestSuccs[i], currentSucc)
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
            if (rand_r(&thSpecific->rndState) > graspThreshold)
                break;
    }

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        bestSuccs[successorSubIndex].cost = computeEdgeCost(thSpecific->X[lastAddedPos], thSpecific->Y[lastAddedPos], thSpecific->X[bestSuccs[successorSubIndex].node], thSpecific->Y[bestSuccs[successorSubIndex].node], inst);
    #endif

    return bestSuccs[successorSubIndex];
}
#endif
