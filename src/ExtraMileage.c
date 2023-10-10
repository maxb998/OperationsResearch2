#include "Tsp.h"


#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

#define AVX_SHIFT 0

typedef struct
{
    Solution bestSol;

    enum EMInitType startOption;

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
static ThreadSharedData initThreadSharedData (Instance *inst, enum EMInitType startOption, double timeLimit);

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
static bool checkSolutionIntegrity(ThreadSpecificData *thSpecific);

// Add node specified in succ into the tour(from 0 to nCovered) inside thSpecific.workingSol without compromising integrity of the solution
static inline void insertNodeInSolution(ThreadSpecificData *thSpecific, int nCovered, SuccessorData succ);

// Find best or close to best successor
static SuccessorData findSuccessor(ThreadSpecificData *thSpecific, int nCovered);


Solution ExtraMileage(Instance *inst, enum EMInitType startOption, double timeLimit, int nThreads)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError("ExtraMileage: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = inst->params.nThreads;
    
    if ((inst->params.graspType == GRASP_NONE) && (startOption == EM_INIT_FARTHEST_POINTS)) // only one solution exists with these settings so only need one thread working
        nThreads = 1;

    // apply extra mileage
    ThreadSharedData thShared = initThreadSharedData(inst, startOption, startTime + timeLimit);

    ThreadSpecificData thSpecifics[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, (unsigned int)rand());
        pthread_create(&threads[i], NULL, runExtraMileage, &thSpecifics[i]);
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

void applyExtraMileage(Solution *sol, int nCovered, unsigned int *rndState)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    // check input solution
    for (int i = 0; i < nCovered; i++)
        for (int j = 0; j < nCovered; j++)
            if ((i != j) && (sol->indexPath[i] == sol->indexPath[j]))
                throwError("applyExtraMileage: Presented solution has the same node in the path before nCovered. sol->indexPath[%d] = sol->indexPath[%d] = %d", i, j, sol->indexPath[i]);

    // recompute cost of nCovered loop
    sol->cost = 0;
    for (int i = 0; i < nCovered-1; i++)
        sol->cost += cvtFloat2Cost(computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i+1]], inst->Y[sol->indexPath[i+1]], inst->params.edgeWeightType, inst->params.roundWeights));
    sol->cost += cvtFloat2Cost(computeEdgeCost(inst->X[sol->indexPath[nCovered-1]], inst->Y[sol->indexPath[nCovered-1]], inst->X[sol->indexPath[0]], inst->Y[sol->indexPath[0]], inst->params.edgeWeightType, inst->params.roundWeights));
            
    // make thSpecific.workingSol consistent with the Extra Mileage implementation(only nodes after nCovered)
    bool *isNodeInSol = calloc(n, sizeof(bool));
    for (int i = 0; i < nCovered; i++)
        isNodeInSol[sol->indexPath[i]] = true;
    for (int i = 0, j = nCovered-1; i < n; i++)
        if (!isNodeInSol[i])
            sol->indexPath[j++] = i;
    free(isNodeInSol);
    for (int i = n; i < n + AVX_VEC_SIZE; i++)
        sol->indexPath[i] = i;

    ThreadSharedData unallocatedThShared = { .bestSol = { .instance = sol->instance } };
    ThreadSpecificData thSpecific = initThreadSpecificData(&unallocatedThShared, *rndState);

    // must clone solution since ThreadSpecificData will free workingSol.indexPath and we don't want to free sol.indexPath
    cloneSolution(sol, &thSpecific.workingSol);
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        for (int i = 0; i < (n + AVX_VEC_SIZE) * 2; i++)
            thSpecific.X[i] = inst->X[sol->indexPath[i]];
    #endif

    applyExtraMileage_Internal(&thSpecific, nCovered);

    cloneSolution(&thSpecific.workingSol, sol);

    *rndState = thSpecific.rndState;
    destroyThreadSpecificData(&thSpecific);
}


static ThreadSharedData initThreadSharedData (Instance *inst, enum EMInitType startOption, double timeLimit)
{
    ThreadSharedData thShared = {
        .startOption = startOption,
        .timeLimit = timeLimit
    };

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
        thSpecific.X = malloc((inst->nNodes + AVX_VEC_SIZE) * 3 * sizeof(int));
        if (thSpecific.X == NULL)
            throwError("ExtraMileage -> initThreadSpecificData: Failed to allocate memory");
        thSpecific.Y = &thSpecific.X[inst->nNodes + AVX_VEC_SIZE];
        thSpecific.costCache = &thSpecific.Y[inst->nNodes + AVX_VEC_SIZE];
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

        applyExtraMileage_Internal(thSpecific, 2);

        if (thSpecific->workingSol.cost < thShared->bestSol.cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol.cost)
                updateBestSolution(thSpecific);
            pthread_mutex_unlock(&thShared->mutex);
        }

        thSpecific->iterCount++;

        if ((inst->params.graspType == GRASP_NONE) && (thShared->startOption == EM_INIT_FARTHEST_POINTS)) // if true only this solution will be found so we can exit
            break;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }

    return NULL;
}

static void initialization(ThreadSpecificData *thSpecific)
{
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;
    int *indexPath = thSpecific->workingSol.indexPath;

    // Set data for solution (copy coords from distance and create index path as 0,1,2,...,n-1)
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        for (int i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
            thSpecific->X[i] = inst->X[i];
    #endif
    for (int i = 0; i < inst->nNodes + AVX_VEC_SIZE; i++)
        indexPath[i] = i;
    // set all costCache to 0 to avoid issues with selection using avx
    for (int i = 0; i < inst->nNodes + AVX_VEC_SIZE; i++)
        thSpecific->costCache[i] = -INFINITY;
    

    switch (thShared->startOption)
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
    if (bestSol->instance->params.logLevel >= LOG_LVL_DEBUG)
        if (!checkSolution(newBest))
		    throwError("updateBestSolutionEM: newBest Solution is not valid");
    
    LOG(LOG_LVL_LOG, "Found better solution: cost = %lf", cvtCost2Double(newBest->cost));

    register int *temp;
    swapElems(bestSol->indexPath, newBest->indexPath, temp);
}

static void applyExtraMileage_Internal(ThreadSpecificData *thSpecific, int nCovered)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
        int *indexPath = thSpecific->workingSol.indexPath;
    #endif

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        enum EdgeWeightType ewt = thSpecific->workingSol.instance->params.edgeWeightType;
        bool roundW = thSpecific->workingSol.instance->params.roundWeights;
    #endif

    if (!checkSolutionIntegrity(thSpecific))
        throwError("applyExtraMileage_Internal: thSpecific.workingSol is not consistent");

    // save element to last position & close the tour at index nCovered
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        thSpecific->X[n] = thSpecific->X[nCovered];
        thSpecific->Y[n] = thSpecific->Y[nCovered];
        thSpecific->X[nCovered] = thSpecific->X[0];
        thSpecific->Y[nCovered] = thSpecific->Y[0];
    #endif
    thSpecific->workingSol.indexPath[n] = thSpecific->workingSol.indexPath[nCovered];
    thSpecific->workingSol.indexPath[nCovered] = thSpecific->workingSol.indexPath[0];

    // setup costCache
    for (int i = 0; i < nCovered; i++)
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            thSpecific->costCache[i] = computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[i+1], thSpecific->Y[i+1], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            thSpecific->costCache[i] = computeEdgeCost(inst->X[indexPath[i]], inst->Y[indexPath[i]], inst->X[indexPath[i+1]], inst->Y[indexPath[i+1]], ewt, roundW);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            thSpecific->costCache[i] = inst->edgeCostMat[indexPath[i] * n + indexPath[i+1]];
        #endif

    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    for (; nCovered < n; nCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        SuccessorData succ;

        if ((inst->params.graspType == GRASP_RANDOM) && (graspThreshold > rand_r(&thSpecific->rndState)))
        {
            succ.node = genRandom(&thSpecific->rndState, (nCovered + 1), (n + 1));
            succ.anchor = genRandom(&thSpecific->rndState, 0, nCovered);

            #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
                succ.newCost0 = computeEdgeCost(thSpecific->X[succ.node], thSpecific->Y[succ.node], thSpecific->X[succ.anchor], thSpecific->Y[succ.anchor], ewt, roundW);
                succ.newCost1 = computeEdgeCost(thSpecific->X[succ.node], thSpecific->Y[succ.node], thSpecific->X[succ.anchor+1], thSpecific->Y[succ.anchor+1], ewt, roundW);
            
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                succ.newCost0 = computeEdgeCost(inst->X[indexPath[succ.node]], inst->Y[indexPath[succ.node]], inst->X[indexPath[succ.anchor]], inst->Y[indexPath[succ.anchor]], ewt, roundW);
                succ.newCost1 = computeEdgeCost(inst->X[indexPath[succ.node]], inst->Y[indexPath[succ.node]], inst->X[indexPath[succ.anchor+1]], inst->Y[indexPath[succ.anchor+1]], ewt, roundW);

            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                succ.newCost0 = inst->edgeCostMat[(size_t)indexPath[succ.node] * (size_t)n + (size_t)indexPath[succ.anchor]];
                succ.newCost1 = inst->edgeCostMat[(size_t)indexPath[succ.node] * (size_t)n + (size_t)indexPath[succ.anchor+1]];

            #endif
        }
        else
            succ = findSuccessor(thSpecific, nCovered);

        insertNodeInSolution(thSpecific, nCovered, succ);
    }
}

static bool checkSolutionIntegrity(ThreadSpecificData *thSpecific)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    for (int i = 0; i < n; i++)
    {
        int index = thSpecific->workingSol.indexPath[i];
        if (index < 0 || index > n)
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: workingSol.indexPath[%lu] = %d which is not within the limits", i, index);
            return false;
        }
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            if ((thSpecific->X[i] != inst->X[index]) || (thSpecific->Y[i] != inst->Y[index]))
            {
                LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: Mismatch at index %lu in solution", i);
                return false;
            }
        #endif
    }

    // everything checks out
    return true;
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

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX) && AVX_SHIFT
        // shift elements forward of 1 position iteratively with avx until vector is too big for the amount of elements to shift (do AVX_VEC_SIZE elements per iteration)
        for (i -= AVX_VEC_SIZE; i > succ.anchor; i -= AVX_VEC_SIZE)
        {
            __m256 xVec = _mm256_loadu_ps(&thSpecific->X[i]);
            __m256 yVec = _mm256_loadu_ps(&thSpecific->Y[i]);
            __m256i indexVec = _mm256_loadu_si256((__m256i_u*)&indexPath[i]);
            _mm256_storeu_ps(&thSpecific->X[i + 1], xVec);
            _mm256_storeu_ps(&thSpecific->Y[i + 1], yVec);
            _mm256_storeu_si256((__m256i_u*)&indexPath[i + 1], indexVec);
        }
        i += AVX_VEC_SIZE;
    #endif

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

    // Contains best cost values (extra cost, current cost, new cost 0 and new cost 1)
    __m256 bestExtraCostVec = _mm256_set1_ps(INFINITY);
    __m256 bestAltEdgeCostVec0 = bestExtraCostVec, bestAltEdgeCostVec1 = bestExtraCostVec;
    // Contains the indexes of the nodes from which the best (chosen according to bestMileageVec) one will be added to the solution at the end of the iteration
    __m256i bestNodesVec = _mm256_set1_epi32(-1);
    // Contains the indexes corresponding to the edge that will be removed/ splitted to accomodate the new node
    __m256i bestAnchorsVec = _mm256_set1_epi32(-1);

    for (int node = nCovered+1; node <= n; node++)
    {
        __m256 xNode = _mm256_broadcast_ss(&thSpecific->X[node]), yNode = _mm256_broadcast_ss(&thSpecific->Y[node]);

        // Vectors to keep track of the best node index and it's anchor index
        __m256i currNodeIndex = _mm256_set1_epi32(node);
        __m256i anchorsIndexesVec = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
        __m256i incrementVec = _mm256_set1_epi32(AVX_VEC_SIZE);

        for (int anchors = 0; anchors < nCovered; anchors+=AVX_VEC_SIZE)
        {
            __m256 currEdgeCostVec = _mm256_loadu_ps(&thSpecific->costCache[anchors]);
            
            __m256 altEdgeCostVec0 = computeEdgeCost_VEC(xNode, yNode, _mm256_loadu_ps(&thSpecific->X[anchors]), _mm256_loadu_ps(&thSpecific->Y[anchors]), ewt, roundW);
            __m256 altEdgeCostVec1 = computeEdgeCost_VEC(xNode, yNode, _mm256_loadu_ps(&thSpecific->X[anchors+1]), _mm256_loadu_ps(&thSpecific->Y[anchors+1]), ewt, roundW);
            __m256 currExtraCostVec = _mm256_sub_ps(_mm256_add_ps(altEdgeCostVec0, altEdgeCostVec1), currEdgeCostVec);

            // Compare extra cost and get binary mask
            __m256 cmpMask = _mm256_cmp_ps(currExtraCostVec, bestExtraCostVec, _CMP_LT_OQ);

            // keep only the best
            bestExtraCostVec = _mm256_blendv_ps(bestExtraCostVec, currExtraCostVec, cmpMask);
            bestAltEdgeCostVec0 = _mm256_blendv_ps(bestAltEdgeCostVec0, altEdgeCostVec0, cmpMask);
            bestAltEdgeCostVec1 = _mm256_blendv_ps(bestAltEdgeCostVec1, altEdgeCostVec1, cmpMask);
            bestAnchorsVec = _mm256_blendv_epi8(bestAnchorsVec, anchorsIndexesVec, _mm256_castps_si256(cmpMask));
            bestNodesVec = _mm256_blendv_epi8(bestNodesVec, currNodeIndex, _mm256_castps_si256(cmpMask));
            
            // move to next operation
            anchorsIndexesVec = _mm256_add_epi32(anchorsIndexesVec, incrementVec);
        }
    }
    // at this point we must select a candidate(best or almost-best(grasp))

    // used to store data from any avx registers to memory
    float avxStoreFloat[AVX_VEC_SIZE];
    _mm256_storeu_ps(avxStoreFloat, bestExtraCostVec);

    int sortedArgs[AVX_VEC_SIZE];
    argsort(avxStoreFloat, sortedArgs, AVX_VEC_SIZE);

    int chosenIndex = sortedArgs[0];
    if ((inst->params.graspType == GRASP_ALMOSTBEST) && (n - nCovered > AVX_VEC_SIZE + 1) && !(avxStoreFloat[sortedArgs[1]] == -INFINITY) && (graspThreshold > rand_r(&thSpecific->rndState)))
        for (int i = 1; i < AVX_VEC_SIZE - 1; i++)
            if ((avxStoreFloat[sortedArgs[i+1]] == -INFINITY) || (rand_r(&thSpecific->rndState) < graspThreshold))
                break;

    SuccessorData retVal;
    int *avxStoreInt = (int*)avxStoreFloat;
    
    _mm256_storeu_si256((__m256i_u *)avxStoreInt, bestNodesVec);
    retVal.node = avxStoreInt[chosenIndex];
    _mm256_storeu_si256((__m256i_u *)avxStoreInt, bestAnchorsVec);
    retVal.anchor = avxStoreInt[chosenIndex];
    _mm256_storeu_ps(avxStoreFloat, bestAltEdgeCostVec0);
    retVal.newCost0 = avxStoreFloat[chosenIndex];
    _mm256_storeu_ps(avxStoreFloat, bestAltEdgeCostVec1);
    retVal.newCost1 = avxStoreFloat[chosenIndex];

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

