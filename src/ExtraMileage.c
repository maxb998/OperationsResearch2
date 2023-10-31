#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

//#define DEBUG

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

    int initIndexes[2];

    float *X;
    float *Y;
    float *costCache;
    Solution workingSol;

} ThreadSpecificData;

typedef struct 
{
    int node;
    int anchor;

    float newCost0;
    float newCost1;

    float extraCost;
} SuccessorData;


// Setup internal variables and initializes mutex
static ThreadSharedData initThreadSharedData (Instance *inst, double timeLimit);

// Destroy mutex
static void destroyThreadSharedData (ThreadSharedData *thShared);

// Setup internal variables and allocate memory in thSpecific.X/Y if using COMPUTE_OPTION_AVX
static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState);

// Deallocate memory pointed by thSpecific.X/Y if necessary
static void destroyThreadSpecificData(ThreadSpecificData *thSpecific);

// Swap elements in thSpecific.workingSol.indexPath and hSpecific.X, thSpecific.Y if COMPUTE_OPTION_AVX
static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2);

// Method called by possibly multiple threads to run extra mileage until time limit. arg is a pointer to a ThreadSpecificData initialized struct
static void *runExtraMileage(void *arg);

// Applies extra mileage to correctly initialized solution(that is thSpecific.workingSol)
static void applyExtraMileage_Internal(ThreadSpecificData *thSpecific, int nCovered);

// Initializes a solution in either a randomized way or by selecting the farthest points
static void initialization(ThreadSpecificData *thSpecific);

// Swap structures thSpecific.thShared.bestSol with thSpecific.workingSol and prints a notification message
static void updateBestSolution(ThreadSpecificData *thSpecific);

// Called by initializations, does exactly what its name says
static void farthestPointsInit(ThreadSpecificData *thSpecific);

// Check whether solution is correctly structured to run extra mileage(after initialization or when initialized solution comes from outside)
static void checkSolutionIntegrity(ThreadSpecificData *thSpecific, int nCovered);

// Add node specified in succ into the tour(from 0 to nCovered) inside thSpecific.workingSol without compromising integrity of the solution
static inline void insertNodeInSolution(ThreadSpecificData *thSpecific, int nCovered, SuccessorData succ);

// Find best or close to best successor
static SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int nCovered);


Solution ExtraMileage(Instance *inst, double timeLimit)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // apply extra mileage
    ThreadSharedData thShared = initThreadSharedData(inst, startTime + timeLimit);

    ThreadSpecificData thSpecifics[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, (unsigned int)rand());
        pthread_create(&threads[i], NULL, runExtraMileage, &thSpecifics[i]);
    }

    int iterCount = 0;
    for (int i = 0; i < inst->params.nThreads; i++)
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

static ThreadSharedData initThreadSharedData (Instance *inst, double timeLimit)
{
    ThreadSharedData thShared = { .timeLimit = timeLimit };

    if (pthread_mutex_init(&thShared.mutex, NULL)) throwError("ExtraMileage -> initThreadSharedData: Failed to initialize mutex");

    thShared.bestSol = newSolution(inst);

    return thShared;
}

static void destroyThreadSharedData (ThreadSharedData *thShared)
{
    if (pthread_mutex_init(&thShared->mutex, NULL)) throwError("ExtraMileage -> destroyThreadSharedData: Failed to destroy mutex");
    destroySolution(&thShared->bestSol);
}

static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState)
{
    ThreadSpecificData thSpecific = {
        .thShared=thShared,
        .rndState=rndState,
        .iterCount=0,
        .initIndexes={-1,-1}
    };

    Instance *inst = thShared->bestSol.instance;
    thSpecific.workingSol=newSolution(inst);

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        thSpecific.X = malloc((inst->nNodes + AVX_VEC_SIZE) * (2 * sizeof(int) + sizeof(float)));
        if (thSpecific.X == NULL)
            throwError("ExtraMileage -> initThreadSpecificData: Failed to allocate memory");
        thSpecific.Y = &thSpecific.X[inst->nNodes + AVX_VEC_SIZE];
        thSpecific.costCache = (float*)&thSpecific.Y[inst->nNodes + AVX_VEC_SIZE];
    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE))
        thSpecific.costCache = malloc((inst->nNodes + AVX_VEC_SIZE) * sizeof(float));
        if (thSpecific.costCache == NULL)
            throwError("ExtraMileage -> initThreadSpecificData: Failed to allocate memory");
    #endif

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        free(thSpecific->X);
    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE))
        free(thSpecific->costCache);
    #endif

    thSpecific->X = thSpecific->Y = thSpecific->costCache = NULL;
}

static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2)
{
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        register float tempFloat;
        swapElems(thSpecific->X[pos1], thSpecific->X[pos2], tempFloat);
        swapElems(thSpecific->Y[pos1], thSpecific->Y[pos2], tempFloat);
    #endif
    register int tempInt;
    swapElems(thSpecific->workingSol.indexPath[pos1], thSpecific->workingSol.indexPath[pos2], tempInt);
}

static void *runExtraMileage(void * arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    while (currentTime < thShared->timeLimit)
    {
        initialization(thSpecific);

        #if ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
            int *indexPath = thSpecific->workingSol.indexPath;
        #endif

        #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            enum EdgeWeightType ewt = thSpecific->workingSol.instance->params.edgeWeightType;
            bool roundW = thSpecific->workingSol.instance->params.roundWeights;
        #endif

        // setup costCache
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            thSpecific->costCache[0] = thSpecific->costCache[1] = computeEdgeCost(thSpecific->X[0], thSpecific->Y[0], thSpecific->X[1], thSpecific->Y[1], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            thSpecific->costCache[0] = thSpecific->costCache[1] = computeEdgeCost(inst->X[indexPath[0]], inst->Y[indexPath[0]], inst->X[indexPath[1]], inst->Y[indexPath[1]], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            thSpecific->costCache[0] = thSpecific->costCache[1] = inst->edgeCostMat[indexPath[0] * inst->nNodes + indexPath[1]];
        #endif

        applyExtraMileage_Internal(thSpecific, 2);

        if (thSpecific->workingSol.cost < thShared->bestSol.cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol.cost)
                updateBestSolution(thSpecific);
            pthread_mutex_unlock(&thShared->mutex);
        }

        thSpecific->iterCount++;

        if ((inst->params.graspType == GRASP_NONE) && (inst->params.nThreads == EM_INIT_FARTHEST_POINTS)) // if true only this solution will be found so we can exit
            break;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }

    return NULL;
}

static void initialization(ThreadSpecificData *thSpecific)
{
    Instance *inst = thSpecific->workingSol.instance;
    int *indexPath = thSpecific->workingSol.indexPath;

    // Set data for solution (copy coords from distance and create index path as 0,1,2,...,n-1)
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        for (int i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
            thSpecific->X[i] = inst->X[i];
    #endif
    for (int i = 0; i < inst->nNodes + AVX_VEC_SIZE; i++)
        indexPath[i] = i;
    #if (COMPUTATION_TYPE != COMPUTE_OPTION_AVX)
        // set all costCache to 0 to avoid issues with selection using avx
        for (int i = 0; i < inst->nNodes + AVX_VEC_SIZE; i++)
            thSpecific->costCache[i] = -INFINITY;
    #endif
    

    switch (inst->params.emInitOption)
    {
    case EM_INIT_RANDOM:
        {
        // select two random nodes
        thSpecific->initIndexes[0] = genRandom(&thSpecific->rndState, 0, inst->nNodes), thSpecific->initIndexes[1] = genRandom(&thSpecific->rndState, 0, inst->nNodes);
        while (thSpecific->initIndexes[0] == thSpecific->initIndexes[1])
            thSpecific->initIndexes[1] = genRandom(&thSpecific->rndState, 0, inst->nNodes);
        }
        break;

    case EM_INIT_FARTHEST_POINTS:
        if (thSpecific->initIndexes[0] == -1)
            farthestPointsInit(thSpecific);
        break;
    }

    // add the two nodes to the solution (order does not matter)
    swapElemsInThSpecific(thSpecific, 0, thSpecific->initIndexes[0]);
    swapElemsInThSpecific(thSpecific, 1, thSpecific->initIndexes[1]);

    // update cost
    float firstEdgeCost;
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        firstEdgeCost = computeEdgeCost(inst->X[indexPath[0]], inst->Y[indexPath[0]], inst->X[indexPath[1]], inst->Y[indexPath[1]], inst->params.edgeWeightType , inst->params.roundWeights);
        thSpecific->workingSol.cost = cvtFloat2Cost(firstEdgeCost) * 2;
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        firstEdgeCost = inst->edgeCostMat[(size_t)indexPath[0] * (size_t)inst->nNodes + (size_t)indexPath[1]];
        thSpecific->workingSol.cost = cvtFloat2Cost(firstEdgeCost) * 2;
    #endif
}

static inline void farthestPointsInit(ThreadSpecificData *thSpecific)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        __m256 maxCostVec = _mm256_set1_ps(0), rowMaxCostVec = _mm256_set1_ps(0); // cost is always positive
        __m256i maxIndexVec1 = _mm256_set1_epi32(0), maxIndexVec2 = _mm256_set1_epi32(0);
        __m256i incrementVec = _mm256_set1_epi32(AVX_VEC_SIZE), ones = _mm256_set1_epi32(1);

        __m256i rowIDsVec = _mm256_set1_epi32(0); // the content of this are always all i
        for (int i = 0; i < n - 1; i++, rowIDsVec = _mm256_add_epi32(rowIDsVec, ones))
        {
            __m256i colIDsVec = _mm256_set_epi32(8 + i, 7 + i, 6 + i, 5 + i, 4 + i, 3 + i, 2 + i, 1 + i);
            __m256 x1 = _mm256_broadcast_ss(&inst->X[i]), y1 = _mm256_broadcast_ss(&inst->Y[i]);

            for (int j = i + 1; j < n; j += AVX_VEC_SIZE, colIDsVec = _mm256_add_epi32(colIDsVec, incrementVec))
            {
                if (j > n - AVX_VEC_SIZE)
                {
                    colIDsVec = _mm256_sub_epi32(colIDsVec, _mm256_set1_epi32(AVX_VEC_SIZE - n % AVX_VEC_SIZE));
                    j = n - AVX_VEC_SIZE; // subtract AVX_VEC_SIZE one extra time to compensate the increment of the loop
                }

                __m256 x2 = _mm256_loadu_ps(&inst->X[j]), y2 = _mm256_loadu_ps(&inst->Y[j]);
                __m256 costVec = computeEdgeCost_VEC(x1, y1, x2, y2, inst->params.edgeWeightType, inst->params.roundWeights);

                // check if there are costier connections in this iteration and save results
                __m256 mask = _mm256_cmp_ps(costVec, rowMaxCostVec, _CMP_GT_OQ);
                rowMaxCostVec = _mm256_blendv_ps(rowMaxCostVec, costVec, mask);
                maxIndexVec1 = _mm256_blendv_epi8(maxIndexVec1, colIDsVec, _mm256_castps_si256(mask));
            }

            __m256 mask = _mm256_cmp_ps(rowMaxCostVec, maxCostVec, _CMP_GT_OQ);
            maxCostVec = _mm256_blendv_ps(maxCostVec, rowMaxCostVec, mask);
            maxIndexVec2 = _mm256_blendv_epi8(maxIndexVec2, rowIDsVec, _mm256_castps_si256(mask));
        }

        float maxCosts[AVX_VEC_SIZE];
        _mm256_storeu_ps(maxCosts, maxCostVec);
        int maxIndex = 0;
        for (int i = 1; i < AVX_VEC_SIZE; i++)
            if (maxCosts[maxIndex] < maxCosts[i])
                maxIndex = i;

        int maxIndexes[AVX_VEC_SIZE];
        _mm256_storeu_si256((__m256i_u *)maxIndexes, maxIndexVec1);
        thSpecific->initIndexes[0] = maxIndexes[maxIndex];
        _mm256_storeu_si256((__m256i_u *)maxIndexes, maxIndexVec2);
        thSpecific->initIndexes[1] = maxIndexes[maxIndex];
    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
        float maxCost = 0;
        for (int i = 0; i < n - 1; i++)
        {
            #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                // hope the compiler loads into the registers these variables that are accessed every time in the loop
                float x1 = inst->X[i];
                float y1 = inst->Y[i];
            #endif

            for (int j = i + 1; j < n; j++)
            {
                float currentCost;
                #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                    currentCost = computeEdgeCost(x1, y1, inst->X[j], inst->Y[j], inst->params.edgeWeightType, inst->params.roundWeights);
                #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                    currentCost = inst->edgeCostMat[(size_t)i * (size_t)n + (size_t)j];
                #endif

                if (currentCost > maxCost)
                {
                    maxCost = currentCost;
                    thSpecific->initIndexes[0] = i;
                    thSpecific->initIndexes[1] = j;
                }
            }
        }
    #endif

    LOG(LOG_LVL_DEBUG, "Extra Mileage EM_INIT_FARTHEST_POINTS: maximum cost found between nodes %d and %d", thSpecific->initIndexes[0], thSpecific->initIndexes[1]);
}

static void updateBestSolution(ThreadSpecificData *thSpecific)
{
    Solution *bestSol = &thSpecific->thShared->bestSol;
    Solution *newBest = &thSpecific->workingSol;

    // swap cost before anything else to avoid making threads wait mutex more than they should
    bestSol->cost = newBest->cost;

    // check solution when debugging
    #ifdef DEBUG
        if (!checkSolution(newBest))
		    throwError("updateBestSolutionEM: newBest Solution is not valid");
    #endif
    
    LOG(LOG_LVL_LOG, "Found better solution: cost = %lf", cvtCost2Double(newBest->cost));

    register int *temp;
    swapElems(bestSol->indexPath, newBest->indexPath, temp);
}

static void applyExtraMileage_Internal(ThreadSpecificData *thSpecific, int nCovered)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        enum EdgeWeightType ewt = thSpecific->workingSol.instance->params.edgeWeightType;
        bool roundW = thSpecific->workingSol.instance->params.roundWeights;
    #endif

    // save element to last position & close the tour at index nCovered
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        thSpecific->X[n] = thSpecific->X[nCovered];
        thSpecific->Y[n] = thSpecific->Y[nCovered];
        thSpecific->X[nCovered] = thSpecific->X[0];
        thSpecific->Y[nCovered] = thSpecific->Y[0];
    #endif
    thSpecific->workingSol.indexPath[n] = thSpecific->workingSol.indexPath[nCovered];
    thSpecific->workingSol.indexPath[nCovered] = thSpecific->workingSol.indexPath[0];

    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    for (; nCovered < n; nCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        #ifdef DEBUG
            checkSolutionIntegrity(thSpecific, nCovered);
        #endif

        SuccessorData succ;

        if ((inst->params.graspType == GRASP_RANDOM) && (graspThreshold > rand_r(&thSpecific->rndState)))
        {
            succ.node = genRandom(&thSpecific->rndState, (nCovered + 1), (n + 1));
            succ.anchor = genRandom(&thSpecific->rndState, 0, nCovered);

            #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
                succ.newCost0 = computeEdgeCost(thSpecific->X[succ.node], thSpecific->Y[succ.node], thSpecific->X[succ.anchor], thSpecific->Y[succ.anchor], ewt, roundW);
                succ.newCost1 = computeEdgeCost(thSpecific->X[succ.node], thSpecific->Y[succ.node], thSpecific->X[succ.anchor+1], thSpecific->Y[succ.anchor+1], ewt, roundW);
            
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                int *indexPath = thSpecific->workingSol.indexPath;
                succ.newCost0 = computeEdgeCost(inst->X[indexPath[succ.node]], inst->Y[indexPath[succ.node]], inst->X[indexPath[succ.anchor]], inst->Y[indexPath[succ.anchor]], ewt, roundW);
                succ.newCost1 = computeEdgeCost(inst->X[indexPath[succ.node]], inst->Y[indexPath[succ.node]], inst->X[indexPath[succ.anchor+1]], inst->Y[indexPath[succ.anchor+1]], ewt, roundW);

            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                int *indexPath = thSpecific->workingSol.indexPath;
                succ.newCost0 = inst->edgeCostMat[(size_t)indexPath[succ.node] * (size_t)n + (size_t)indexPath[succ.anchor]];
                succ.newCost1 = inst->edgeCostMat[(size_t)indexPath[succ.node] * (size_t)n + (size_t)indexPath[succ.anchor+1]];

            #endif
        }
        else
            succ = findSuccessor(thSpecific, nCovered);

        insertNodeInSolution(thSpecific, nCovered, succ);
    }
}

static void checkSolutionIntegrity(ThreadSpecificData *thSpecific, int nCovered)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    bool *found = malloc(n * sizeof(bool));
    for (int i = 0; i < n; i++)
        found[i] = false;
    

    for (int i = 0; i < nCovered; i++)
    {
        if (found[thSpecific->workingSol.indexPath[i]])
            throwError("ExtraMileage-checkSolutionIntegrity: Solution is incoherent in the covered part. Node %d appear twice up to position %d with %d covered nodes", thSpecific->workingSol.indexPath[i], i, nCovered);
        found[thSpecific->workingSol.indexPath[i]] = true;
    }

    if (thSpecific->workingSol.indexPath[0] != thSpecific->workingSol.indexPath[nCovered])
        throwError("ExtraMileage-checkSolutionIntegrity: Solution path is incoherent, last position does not match the first one");

    for (int i = nCovered+1; i <= n; i++)
    {
        if (found[thSpecific->workingSol.indexPath[i]])
            throwError("ExtraMileage, checkSolutionIntegrity: Solution is incoherent. Node %d appear twice up to position %d with %d covered nodes", thSpecific->workingSol.indexPath[i], i, nCovered);
        found[thSpecific->workingSol.indexPath[i]] = true;
    }

    // check if all nodes are present
    for (int i = 0; i < n; i++)
        if (!found[i])
            throwError("ExtraMileage-checkSolutionIntegrity: Node %d was not present in the solution", i);
    free(found);
    
    // check X and Y
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        for (int i = 0; i <= n; i++)
        {
            if (thSpecific->X[i] != inst->X[thSpecific->workingSol.indexPath[i]])
                throwError("ExtraMileage-checkSolutionIntegrity: thSpecific.X[%d] does not match correctly", i);
            if (thSpecific->Y[i] != inst->Y[thSpecific->workingSol.indexPath[i]])
                throwError("ExtraMileage-checkSolutionIntegrity: thSpecific.Y[%d] does not match correctly", i);
        }
    #endif

    // check cost and costCache
    __uint128_t recomputedCost = 0;
    for (int i = 0; i < nCovered; i++)
    {
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            enum EdgeWeightType ewt = inst->params.edgeWeightType;
            bool roundW = inst->params.roundWeights;
            float currEdgeCost = computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[i+1], thSpecific->Y[i+1], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            enum EdgeWeightType ewt = inst->params.edgeWeightType;
            bool roundW = inst->params.roundWeights;
            int *indexPath = thSpecific->workingSol.indexPath;
            float currEdgeCost = computeEdgeCost(inst->X[indexPath[i]], inst->Y[indexPath[i]], inst->X[indexPath[i+1]], inst->Y[indexPath[i+1]], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            int *indexPath = thSpecific->workingSol.indexPath;
            float currEdgeCost = inst->edgeCostMat[(size_t)indexPath[i] * (size_t)n + (size_t)indexPath[i+1]];
        #endif
        if (currEdgeCost != thSpecific->costCache[i])
            throwError("ExtraMileage-checkSolutionIntegrity: costCache is incoherent at position %d", i);
        recomputedCost += cvtFloat2Cost(currEdgeCost);
    }
    if (thSpecific->workingSol.cost != recomputedCost)
        throwError("ExtraMileage-checkSolutionIntegrity: Cost value is incorrect");
}

static inline void insertNodeInSolution(ThreadSpecificData *thSpecific, int nCovered, SuccessorData succ)
{
    int *indexPath = thSpecific->workingSol.indexPath;

    nCovered++;

    thSpecific->workingSol.cost -= cvtFloat2Cost(thSpecific->costCache[succ.anchor]);
    thSpecific->workingSol.cost += cvtFloat2Cost(succ.newCost0);
    thSpecific->workingSol.cost += cvtFloat2Cost(succ.newCost1);

    // save best value
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        float bestX, bestY;
        bestX = thSpecific->X[succ.node];
        bestY = thSpecific->Y[succ.node];
    #endif
    int bestIndex = indexPath[succ.node];

    // place elements to insert in the tour at the end of the covered nodes "set"
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        thSpecific->X[succ.node] = thSpecific->X[nCovered];
        thSpecific->Y[succ.node] = thSpecific->Y[nCovered];
    #endif
    indexPath[succ.node] = indexPath[nCovered];

    int i = nCovered;

    // shift elements forward one at a time
    for (i--; i > succ.anchor; i--)
    {
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            thSpecific->X[i+1] = thSpecific->X[i];
            thSpecific->Y[i+1] = thSpecific->Y[i];
        #endif
        indexPath[i+1] = indexPath[i];
        thSpecific->costCache[i+1] = thSpecific->costCache[i];
    }

    i++;

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        thSpecific->X[i] = bestX;
        thSpecific->Y[i] = bestY;
    #endif
    indexPath[i] = bestIndex;
    thSpecific->costCache[i-1] = succ.newCost0;
    thSpecific->costCache[i] = succ.newCost1;

    LOG(LOG_LVL_EVERYTHING, "Extra Mileage Solution Update: Node %d added to solution between nodes %d and %d", indexPath[i], indexPath[i-1], indexPath[i+1]);
}

#if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
static SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int nCovered)
{
    // shortcuts/decluttering
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundW = inst->params.roundWeights;
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    // Contains best mileage values
    __m256 bestExtraCostVec = _mm256_set1_ps(INFINITY);
    // Contains the indexes of the nodes from which the best (chosen according to bestMileageVec) one will be added to the solution at the end of the iteration
    __m256i bestNodesVec = _mm256_set1_epi32(-1);
    // Contains the indexes corresponding to the edge that will be removed/ splitted to accomodate the new node
    __m256i bestAnchorsVec = _mm256_set1_epi32(-1);

    // we do this to avoid the need of checking the last elements loaded by _mm256_loadu -> exploit the "INFINITY" placed at the end of the last elements in sol.X and sol.Y
    for (int i = 0; i < nCovered; i++)
    {
        // Create vectors containig necessary data on the points attached to the edge i
        __m256 x1Vec = _mm256_broadcast_ss(&thSpecific->X[i]), y1Vec = _mm256_broadcast_ss(&thSpecific->Y[i]);
        __m256 x2Vec = _mm256_broadcast_ss(&thSpecific->X[i + 1]), y2Vec = _mm256_broadcast_ss(&thSpecific->Y[i + 1]);

        // Vector that contains only the cost of the current edge
        __m256 curEdgeCostVec = _mm256_broadcast_ss(&thSpecific->costCache[i]); //computeEdgeCost_VEC(x1Vec, y1Vec, x2Vec, y2Vec, ewt, roundW);

        // Vector that contains only the index of the current edge
        __m256i curEdgeID = _mm256_set1_epi32(i);

        // Vector that keeps track of the IDs of the best candidates for the current edge
        __m256i idsVec = _mm256_add_epi32(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0), _mm256_set1_epi32(nCovered + 1));
        __m256i incrementVec = _mm256_set1_epi32(AVX_VEC_SIZE);

        // check for each edge which ones are the best
        for (int u = nCovered + 1; u <= n; u += AVX_VEC_SIZE, idsVec = _mm256_add_epi32(idsVec, incrementVec))
        {
            __m256 curExtraMileageVec;
            {
                __m256 xuVec = _mm256_loadu_ps(&thSpecific->X[u]), yuVec = _mm256_loadu_ps(&thSpecific->Y[u]);
                __m256 altEdge1CostVec = computeEdgeCost_VEC(xuVec, yuVec, x1Vec, y1Vec, ewt, roundW);
                __m256 altEdge2CostVec = computeEdgeCost_VEC(xuVec, yuVec, x2Vec, y2Vec, ewt, roundW);
                __m256 altEdgeCostVec = _mm256_add_ps(altEdge1CostVec, altEdge2CostVec);
                curExtraMileageVec = _mm256_sub_ps(altEdgeCostVec, curEdgeCostVec);
            }

            // Compare curExtraMileageCostVec with bestExtraMileageVec
            __m256 cmpMask = _mm256_cmp_ps(curExtraMileageVec, bestExtraCostVec, _CMP_LT_OQ);

            // Set new best according to comparison result
            bestExtraCostVec = _mm256_blendv_ps(bestExtraCostVec, curExtraMileageVec, cmpMask);
            bestAnchorsVec = _mm256_blendv_epi8(bestAnchorsVec, curEdgeID, _mm256_castps_si256(cmpMask));
            bestNodesVec = _mm256_blendv_epi8(bestNodesVec, idsVec, _mm256_castps_si256(cmpMask));
        }
    }
    // at this point we must select a candidate(best or almost-best(grasp))

    // used to store data from any avx registers to memory
    float avxStoreFloat[AVX_VEC_SIZE];
    _mm256_storeu_ps(avxStoreFloat, bestExtraCostVec);

    int sortedArgs[AVX_VEC_SIZE];
    argsort(avxStoreFloat, sortedArgs, AVX_VEC_SIZE);

    int chosenIndex = sortedArgs[0];
    if ((inst->params.graspType == GRASP_ALMOSTBEST) && (n - nCovered > AVX_VEC_SIZE + 1)&& (graspThreshold > rand_r(&thSpecific->rndState)))
        for (int i = 1; i < AVX_VEC_SIZE - 1; i++)
            if ((avxStoreFloat[sortedArgs[i+1]] == -INFINITY) || (rand_r(&thSpecific->rndState) < graspThreshold))
                break;

    SuccessorData retVal;
    int *avxStoreInt = (int*)avxStoreFloat;
    
    _mm256_storeu_si256((__m256i_u *)avxStoreInt, bestNodesVec);
    retVal.node = avxStoreInt[chosenIndex];
    _mm256_storeu_si256((__m256i_u *)avxStoreInt, bestAnchorsVec);
    retVal.anchor = avxStoreInt[chosenIndex];
    retVal.newCost0 = computeEdgeCost(thSpecific->X[retVal.node], thSpecific->Y[retVal.node], thSpecific->X[retVal.anchor], thSpecific->Y[retVal.anchor], ewt, roundW);
    retVal.newCost1 = computeEdgeCost(thSpecific->X[retVal.node], thSpecific->Y[retVal.node], thSpecific->X[retVal.anchor+1], thSpecific->Y[retVal.anchor+1], ewt, roundW);

    return retVal;
}

#elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
static SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int nCovered)
{
    // shortcuts/decluttering
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    int *indexPath = thSpecific->workingSol.indexPath;
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        enum EdgeWeightType ewt = inst->params.edgeWeightType ;
        bool roundW = inst->params.roundWeights;
    #endif

    SuccessorData bestSuccs[BASE_GRASP_BEST_SAVE_BUFFER_SIZE];
    for (int i = 0; i < BASE_GRASP_BEST_SAVE_BUFFER_SIZE; i++)
        bestSuccs[i].extraCost = INFINITY;
    
    
    for (int node = nCovered+1; node <= n; node++)
    {
        float altEdgeCost0, altEdgeCost1;
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            altEdgeCost1 = computeEdgeCost(inst->X[indexPath[0]], inst->Y[indexPath[0]], inst->X[indexPath[node]], inst->Y[indexPath[node]], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            altEdgeCost1 = inst->edgeCostMat[indexPath[node] * n + indexPath [0]]; 
        #endif

        for (int anchor = 0; anchor < nCovered; anchor++)
        {
            altEdgeCost0 = altEdgeCost1;
            #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                altEdgeCost1 = computeEdgeCost(inst->X[indexPath[anchor+1]], inst->Y[indexPath[anchor+1]], inst->X[indexPath[node]], inst->Y[indexPath[node]], ewt, roundW);
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                altEdgeCost1 = inst->edgeCostMat[indexPath[node] * n + indexPath[anchor+1]];
            #endif
            
            SuccessorData currSucc = {
                 .node=node,
                 .anchor=anchor, 
                 .extraCost = altEdgeCost0 + altEdgeCost1 - thSpecific->costCache[anchor],
                 .newCost0 = altEdgeCost0,
                 .newCost1 = altEdgeCost1
            };

            for (int i = 0; i < BASE_GRASP_BEST_SAVE_BUFFER_SIZE; i++)
            {
                if (currSucc.extraCost < bestSuccs[i].extraCost)
                {
                    SuccessorData temp;
                    swapElems(currSucc, bestSuccs[i], temp);
                }
            }
        }   
    }
    // reached this point we found the absolute #BASE_GRASP_BEST_SAVE_BUFFER_SIZE best combination of anchors and uncovered nodes possible for this iteration
    int chosedNodeSubIndex = 0;
    if ((inst->params.graspType == GRASP_ALMOSTBEST) && (graspThreshold > rand_r(&thSpecific->rndState)) && (n - nCovered - 1 > BASE_GRASP_BEST_SAVE_BUFFER_SIZE))
    {
        chosedNodeSubIndex = 1;
        for (; chosedNodeSubIndex < BASE_GRASP_BEST_SAVE_BUFFER_SIZE - 1; chosedNodeSubIndex++)
            if (rand_r(&thSpecific->rndState) > RAND_MAX / 2)
                break;
    }

    return bestSuccs[chosedNodeSubIndex];
}
#endif

