#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

#define DEBUG

// Defines how big the tenure should be compared to the number of nodes
#define TENURE_RATIO 1/3
// Defines how many non-improving iterations tabu must do before restarting from the best sol
#define TABU__RESTART_FROM_BEST__THRESHOLD 10000


typedef struct
{
    Solution *bestSol;

    pthread_mutex_t mutex;
    double timeLimit;

    int tenureSize;

} ThreadSharedData;

typedef struct
{
    ThreadSharedData *thShared;
    unsigned int rndState;
    int iterCount;

    Solution workingSol;
    float *X;
    float *Y;
    float *costCache;
    
    int *tenure;
    // keeps the cost of the edge in the tenure backed up, not to improve performance, but to check correctness
    float *costBackup;

    int nextTenurePos;

} ThreadSpecificData;

// Initializes mutex and variables
static ThreadSharedData initThreadSharedData (Solution *sol, double timeLimit);
// Destroy mutex
static void destroyThreadSharedData (ThreadSharedData *thShared);
// Initializes variables (and allocates memory on .X, .Y when using COMPUTE_OPTION_AVX)
static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState);
// Frees memory allocated in X and Y if any
static void destroyThreadSpecificData(ThreadSpecificData *thSpecific);
// Method each thread run when executing Tabu Search until time limit
static void *runTabu(void *arg);
// Set the workingSol equal to thShared.bestSol and clears the tabu tenure(and rearranges .X, .Y accordingly when using COMPUTE_OPTION_AVX)
static void setupThSpecificOnBestSol(ThreadSpecificData *thSpecific);
// Select an edge different than the "forbiddenNode" and which is not in the tenure in a random way.
static inline int randomlySelectEdgeOutsideTenure(ThreadSpecificData *thSpecific, int forbiddenNode);
// Add node thSpecific.workingSol.indexPath[indexPathIndex] to the tenure, which means also removing one when the tenure is full.
static void addNodeToTenure(ThreadSpecificData *thSpecific, int indexPathIndex);
// Perform non-improving 2opt move on thSpecific data using node0 and node1
static inline void performNonImproving2OptMove(ThreadSpecificData *thSpecific, int edge0, int edge1);


void TabuSearch(Solution *sol, double timeLimit, int nThreads)
{
    Instance *inst = sol->instance;

    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // must save solution time since 2opt is gonna increase it inconsistently with vns
    double initialSolutionRuntime = sol->execTime;

    // reset seed if debugging
    if (inst->params.logLevel == LOG_LVL_DEBUG)
        srand(inst->params.randomSeed);

    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError("VariableNeighborhood: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = inst->params.nThreads;

    if (!checkSolution(sol))
        throwError("VariableNeighborhood: Input solution is not valid");

    sol->indexPath[inst->nNodes] = sol->indexPath[0];

    apply2OptBestFix(sol); // isn't necessary since it's solution should already be 2-optimized, but just to be sure

    ThreadSharedData thShared = initThreadSharedData(sol, startTime + timeLimit);
    ThreadSpecificData thSpecifics[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, (unsigned int)rand());
        pthread_create(&threads[i], NULL, runTabu, &thSpecifics[i]);
    }

    int iterCount = 0;
    for (int i = 0; i < nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        iterCount += thSpecifics[i].iterCount;
        destroyThreadSpecificData(&thSpecifics[i]);
    }

    destroyThreadSharedData(&thShared);   

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double execTime = cvtTimespec2Double(timeStruct) - startTime;
    sol->execTime = execTime + initialSolutionRuntime;

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %d", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/execTime);
}

static ThreadSharedData initThreadSharedData (Solution *sol, double timeLimit)
{
    ThreadSharedData thShared = { .timeLimit = timeLimit, .bestSol=sol, .tenureSize= sol->instance->nNodes * TENURE_RATIO };

    if (pthread_mutex_init(&thShared.mutex, NULL)) throwError("VariableNeighborhoodSearch -> initThreadSharedData: Failed to initialize mutex");

    return thShared;
}

static void destroyThreadSharedData (ThreadSharedData *thShared)
{
    if (pthread_mutex_init(&thShared->mutex, NULL)) throwError("VariableNeighborhoodSearch -> destroyThreadSharedData: Failed to destroy mutex");
}

static ThreadSpecificData initThreadSpecificData (ThreadSharedData *thShared, unsigned int rndState)
{
    Instance *inst = thShared->bestSol->instance;
    int n = inst->nNodes;

    ThreadSpecificData thSpecific = {
        .thShared=thShared,
        .rndState=rndState,
        .iterCount=0,
        .nextTenurePos=0
    };

    thSpecific.tenure = malloc(thShared->tenureSize * sizeof(int) * 2);
    if (!thSpecific.tenure)
        throwError("Tabu -> initThreadSpecificData: Failed to allocate memory");
    thSpecific.costBackup = (float*)&thSpecific.tenure[thShared->tenureSize];
    // set all elements in tenure to -1 (empty tenure)
    for (int i = 0; i < thShared->tenureSize; i++)
        thSpecific.tenure[i] = -1;

    thSpecific.workingSol=newSolution(inst);
    
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        thSpecific.X = malloc((n + AVX_VEC_SIZE) * 3 * sizeof(int));
        if (!thSpecific.X)
            throwError("Tabu -> initThreadSpecificData: Failed to allocate memory");
        thSpecific.Y = &thSpecific.X[n + AVX_VEC_SIZE];
        thSpecific.costCache = &thSpecific.Y[n + AVX_VEC_SIZE];
    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
        thSpecific.X = thSpecific.Y = NULL;
        thSpecific.costCache = malloc((n + AVX_VEC_SIZE) * sizeof(float));
        if (!thSpecific.costCache)
            throwError("Tabu -> initThreadSpecificData: Failed to allocate memory");
    #endif

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    free(thSpecific->tenure);
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        free(thSpecific->X);
    #else
        free(thSpecific->costCache);
    #endif

    thSpecific->tenure = NULL;
    thSpecific->X = thSpecific->Y = thSpecific->costCache = NULL;
}

static void *runTabu(void *arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    int nonImprovingIterCount = 0;

    while (currentTime < thShared->timeLimit)
    {
        // use 2opt to optimize (setting edges in the costCache to -INFINITY effectively lock that edges)
        apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->X, thSpecific->Y, thSpecific->costCache);

        if (nonImprovingIterCount > TABU__RESTART_FROM_BEST__THRESHOLD)
        {
            setupThSpecificOnBestSol(thSpecific);
            nonImprovingIterCount = 0;
        }
        if (thSpecific->workingSol.cost < thShared->bestSol->cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol->cost)
            {
                LOG(LOG_LVL_LOG, "Found better solution: cost = %lf", cvtCost2Double(thSpecific->workingSol.cost));
                cloneSolution(&thSpecific->workingSol, thShared->bestSol);
            }
            else
                nonImprovingIterCount++;
            pthread_mutex_unlock(&thShared->mutex);
        }
        else
            nonImprovingIterCount++;

        // now perform a non optimal 2opt move, lock one of the edges involved in the move and add it to the tenure
        int edge0 = randomlySelectEdgeOutsideTenure(thSpecific, -1);
        int edge1 = randomlySelectEdgeOutsideTenure(thSpecific, edge0);

        performNonImproving2OptMove(thSpecific, edge0, edge1);

        int edgeToLock = edge0;
        if (genRandom(&thSpecific->rndState, 0, 10) < 5)
            edgeToLock = edge1;
        
        addNodeToTenure(thSpecific, edgeToLock);

        thSpecific->iterCount++;
        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }

    return NULL;
}

static void setupThSpecificOnBestSol(ThreadSpecificData *thSpecific)
{
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    // clone of best sol must be done with locked mutex otherwise there may be synchronization errors when starting a lot of threads
    pthread_mutex_lock(&thShared->mutex);
    cloneSolution(thShared->bestSol, &thSpecific->workingSol);
    pthread_mutex_unlock(&thShared->mutex);

    #if (COMPUTATION_TYPE == COMPUTATE_OPTION_AVX)
        for (int i = 0; i < n; i++)
        {
            thSpecific->X[i] = inst->X[thSpecific->workingSol.indexPath[i]];
            thSpecific->Y[i] = inst->Y[thSpecific->workingSol.indexPath[i]];
        }
        thSpecific->X[n] = thSpecific->X[0];
        thSpecific->Y[n] = thSpecific->Y[0];

        for (int i = n+1; i < n + AVX_VEC_SIZE; i++)
        {
            thSpecific->X[i] = INFINITY;
            thSpecific->Y[i] = INFINITY;
        }
    #endif

    int *path = thShared->bestSol->indexPath;
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        for (int i = 0; i < n; i++) // build cost cache
            thSpecific->costCache[i] = computeEdgeCost(inst->X[path[i]], inst->Y[path[i]], inst->X[path[i+1]], inst->Y[path[i+1]],
                                                        inst->params.edgeWeightType, inst->params.roundWeights);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        for (int i = 0; i < n; i++) // build cost cache
            thSpecific->costCache[i] = inst->edgeCostMat[path[i] * n + path[i+1]];
    #endif

    for (int i = n+1; i < n + AVX_VEC_SIZE; i++)
        thSpecific->costCache[i] = INFINITY;
    
    // reset tenure
    for (int i = 0; i < thShared->tenureSize; i++)
        thSpecific->tenure[i] = -1;
}

static inline int randomlySelectEdgeOutsideTenure(ThreadSpecificData *thSpecific, int forbiddenNode)
{
    int edge, n = thSpecific->workingSol.instance->nNodes;
    // cannot choose any edge in the tenure and edge with the same value as 
    bool rndValueNotValid = true;
    while (rndValueNotValid)
    {
        edge = genRandom(&thSpecific->rndState, 0, n);
        while(edge == forbiddenNode)
            edge = genRandom(&thSpecific->rndState, 0, n);

        int i = 0;
        for (; i < thSpecific->thShared->tenureSize; i++)
            if (edge == thSpecific->tenure[i])
                break;
        if (i == thSpecific->thShared->tenureSize)
            rndValueNotValid = false;
    }
    return edge;
}

static void addNodeToTenure(ThreadSpecificData *thSpecific, int indexPathIndex)
{
    #ifdef DEBUG
    if ((indexPathIndex < 0 ) || (indexPathIndex >= thSpecific->workingSol.instance->nNodes))
        throwError("Tabu -> addNodeToTenure: indexPathIndex is not valid : %d", indexPathIndex);
    // check wheter the edge we want to add to the tenure is already inside the tenure
    for (int i = 0; i < thSpecific->thShared->tenureSize; i++)
        if (thSpecific->tenure[i] == thSpecific->workingSol.indexPath[indexPathIndex])
            throwError("Tabu -> addNodeToTenure: node already inside the tenure indexPath[%d] = %d", indexPathIndex, thSpecific->workingSol.indexPath[indexPathIndex]);
    #endif

    int nextTenurePos = thSpecific->nextTenurePos;

    if (thSpecific->tenure[nextTenurePos] != -1) // tenure full: remove last element from tenure and then add the new one
    {
        int n = thSpecific->workingSol.instance->nNodes;
        for (int i = 0; i < n; i++)
        {
            if (thSpecific->workingSol.indexPath[i] == thSpecific->tenure[nextTenurePos])
            {
                thSpecific->costCache[i] = thSpecific->costBackup[nextTenurePos];
                break;
            }
        }
    }

    // add one element to the tenure
    thSpecific->tenure[nextTenurePos] = indexPathIndex;
    thSpecific->costBackup[nextTenurePos] = thSpecific->costCache[indexPathIndex];
    thSpecific->costCache[indexPathIndex] = -INFINITY;

    thSpecific->nextTenurePos++;
    if (thSpecific->nextTenurePos == thSpecific->thShared->tenureSize)
        thSpecific->nextTenurePos = 0;
}

static inline void performNonImproving2OptMove(ThreadSpecificData *thSpecific, int edge0, int edge1)
{
    Solution *sol = &thSpecific->workingSol;

    if (edge0 > edge1)
    {
        register int temp;
        swapElems(edge0, edge1, temp);
    }

    float altEdge0Cost, altEdge1Cost;

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        enum EdgeWeightType ewt = sol->instance->params.edgeWeightType;
        bool roundW = sol->instance->params.roundWeights;
        altEdge0Cost = computeEdgeCost(thSpecific->X[edge0], thSpecific->Y[edge0], thSpecific->X[edge1], thSpecific->Y[edge1], ewt, roundW);
        altEdge1Cost = computeEdgeCost(thSpecific->X[edge0+1], thSpecific->Y[edge0+1], thSpecific->X[edge1+1], thSpecific->Y[edge1+1], ewt, roundW);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        Instance *inst = sol->instance;
        int *indexPath = sol->indexPath;
        enum EdgeWeightType ewt = inst->params.edgeWeightType;
        bool roundW = inst->params.roundWeights;
        altEdge0Cost = computeEdgeCost(inst->X[indexPath[edge0]], inst->Y[indexPath[edge0]], inst->X[indexPath[edge1]], inst->Y[indexPath[edge1]], ewt, roundW);
        altEdge1Cost = computeEdgeCost(inst->X[indexPath[edge0+1]], inst->Y[indexPath[edge0+1]], inst->X[indexPath[edge1+1]], inst->Y[indexPath[edge1+1]], ewt, roundW);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        Instance *inst = sol->instance;
        int *indexPath = sol->indexPath;
        int n = inst->nNodes;
        altEdge0Cost = inst->edgeCostMat[(size_t)indexPath[edge0] * (size_t)n + (size_t)indexPath[edge1]];
        altEdge1Cost = inst->edgeCostMat[(size_t)indexPath[edge0+1] * (size_t)n + (size_t)indexPath[edge1+1]];
    #endif

    // update cost
    sol->cost += cvtFloat2Cost(altEdge0Cost) + cvtFloat2Cost(altEdge1Cost) - cvtFloat2Cost(thSpecific->costCache[edge0]) - cvtFloat2Cost(thSpecific->costCache[edge1]);

    int smallID = edge0 + 1, bigID = edge1;

    while (smallID < bigID)
    {
        register int tempInt;
        swapElems(sol->indexPath[smallID], sol->indexPath[bigID], tempInt);
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            register float tempFloat;
            swapElems(thSpecific->X[smallID], thSpecific->X[bigID], tempFloat);
            swapElems(thSpecific->Y[smallID], thSpecific->Y[bigID], tempFloat);
        #endif

        smallID++;
        bigID--;
    }

    // update cost cache
    thSpecific->costCache[edge0] = altEdge0Cost;
    thSpecific->costCache[edge1] = altEdge1Cost;

    smallID = edge0 + 1;
    bigID = edge1 - 1;

    while (smallID < bigID)
    {
        register float tempFloat;
        swapElems(thSpecific->costCache[smallID], thSpecific->costCache[bigID], tempFloat);

        smallID++;
        bigID--;
    }
}

