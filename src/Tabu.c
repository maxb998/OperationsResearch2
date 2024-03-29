#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

//#define DEBUG

typedef struct
{
    int node0;
    int node1;
} Edge;

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
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
    float *X;
    float *Y;
    #endif
    float *costCache;
    
    Edge *tenure;
    // keeps the cost of the edge in the tenure backed up, not to improve performance, but to check correctness
    float *costBackup;

    int nextTenurePos;
    int lastTenurePos;

} ThreadSpecificData;


// Initializes mutex and variables
static ThreadSharedData initThreadSharedData (Solution *sol, int tenureSize, double timeLimit);
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
static void addEdgeToTenure(ThreadSpecificData *thSpecific, int indexPathIndex);
// Remove oldest element in the tenure if tenure is not empty
static void removeOldestElementFromTenure(ThreadSpecificData *thSpecific);
// Perform non-improving 2opt move on thSpecific data using node0 and node1
static inline void performNonImproving2OptMove(ThreadSpecificData *thSpecific, int edge0, int edge1);
#ifdef DEBUG
// Checks the Tenure and the costCache for correctness.
static void checkThSpecificData(ThreadSpecificData *thSpecific);
// Prints the element currently in tenure ordered from first to last into strOut
static void printTenure(ThreadSpecificData *thSpecific, char *strOut);
#endif



void TabuSearch(Solution *sol, double timeLimit)
{
    Instance *inst = sol->instance;

    if (inst->params.tabuTenureSize == -1) // auto tenure size
        inst->params.tabuTenureSize = 2;
    else if (inst->params.tabuTenureSize >= sol->instance->nNodes)
        throwError("Specified Tenure size[%d] is bigger or equal to the instance size[%d]", inst->params.tabuTenureSize, sol->instance->nNodes);
    else if (inst->params.tabuTenureSize >= sol->instance->nNodes * 97/100)
        LOG(LOG_LVL_WARN, "Specified Tenure size[%d] is big considering the size of the instance[%d]. Tabu might get stuck and not respect the time limit", inst->params.tabuTenureSize, sol->instance->nNodes);

    LOG(LOG_LVL_NOTICE, "Tabu tenure size is set to %d", inst->params.tabuTenureSize);


    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // must save solution time since 2opt is gonna increase it inconsistently with vns
    double initialSolutionRuntime = sol->execTime;

    // reset seed if debugging
    #ifdef DEBUG
        srand(inst->params.randomSeed);

        if (!checkSolution(sol))
            throwError("Tabu Search: Input solution is not valid");
    #endif

    sol->indexPath[inst->nNodes] = sol->indexPath[0];

    apply2OptBestFix(sol); // isn't necessary since it's solution should already be 2-optimized, but just to be sure

    ThreadSharedData thShared = initThreadSharedData(sol, inst->params.tabuTenureSize, startTime + timeLimit);
    ThreadSpecificData thSpecifics[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        thSpecifics[i] = initThreadSpecificData(&thShared, (unsigned int)rand());
        pthread_create(&threads[i], NULL, runTabu, &thSpecifics[i]);
    }

    int iterCount = 0;
    for (int i = 0; i < inst->params.nThreads; i++)
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

static ThreadSharedData initThreadSharedData (Solution *sol, int tenureSize, double timeLimit)
{
    ThreadSharedData thShared = { .timeLimit = timeLimit, .bestSol=sol, .tenureSize=tenureSize };

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
        .nextTenurePos=0,
        .lastTenurePos=0
    };

    thSpecific.tenure = malloc(thShared->tenureSize * (sizeof(Edge) + sizeof(float)));
    if (!thSpecific.tenure)
        throwError("Tabu -> initThreadSpecificData: Failed to allocate memory");
    thSpecific.costBackup = (float*)&thSpecific.tenure[thShared->tenureSize];
    // set all elements in tenure to -1 (empty tenure)
    for (int i = 0; i < thShared->tenureSize; i++)
        thSpecific.tenure[i].node0 = -1;

    thSpecific.workingSol=newSolution(inst);
    
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        thSpecific.costCache = malloc((n + AVX_VEC_SIZE) * 3 * sizeof(int));
        if (!thSpecific.costCache)
            throwError("Tabu -> initThreadSpecificData: Failed to allocate memory");
        thSpecific.X = &thSpecific.costCache[n + AVX_VEC_SIZE];
        thSpecific.Y = &thSpecific.X[n + AVX_VEC_SIZE];
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
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
    free(thSpecific->costCache);

    thSpecific->tenure = NULL;
    thSpecific->costCache = NULL;
    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        thSpecific->X = thSpecific->Y = NULL;
    #endif
}

static void *runTabu(void *arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    setupThSpecificOnBestSol(thSpecific);

    int nonImprovingIterCount = 0;
    int restartThreshold = thSpecific->workingSol.instance->params.metaRestartThreshold;

    while (currentTime < thShared->timeLimit)
    {
        #ifdef DEBUG
            char tenureStr[3000];
            printTenure(thSpecific, tenureStr);
            LOG(LOG_LVL_DEBUG, tenureStr);        

            LOG(LOG_LVL_DEBUG, "Applied 2otp");
            checkThSpecificData(thSpecific);
        #endif

        // now perform a non optimal 2opt move, lock one of the edges involved in the move and add it to the tenure
        int edge0 = randomlySelectEdgeOutsideTenure(thSpecific, -1);
        int edge1 = randomlySelectEdgeOutsideTenure(thSpecific, edge0);

        performNonImproving2OptMove(thSpecific, edge0, edge1);

        #ifdef DEBUG
            LOG(LOG_LVL_DEBUG, "Performed \"bad\" move");
            checkThSpecificData(thSpecific);
        #endif

        int edgeToLock = edge0;
        if (genRandom(&thSpecific->rndState, 0, 2) == 1)
            edgeToLock = edge1;
        
        addEdgeToTenure(thSpecific, edgeToLock);

        #ifdef DEBUG
            LOG(LOG_LVL_DEBUG, "Added an edge to tenure");
            checkThSpecificData(thSpecific);
        #endif

        int n2OptMoves;
        // use 2opt to optimize (setting edges in the costCache to -INFINITY effectively lock that edges)
        #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            n2OptMoves = apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->X, thSpecific->Y, thSpecific->costCache);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            n2OptMoves = apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->costCache);
        #endif

        while (n2OptMoves > 0)
        {
            removeOldestElementFromTenure(thSpecific);
            #ifdef DEBUG
                LOG(LOG_LVL_DEBUG, "Removed an edge from tenure");
                checkThSpecificData(thSpecific);
            #endif

            #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                n2OptMoves = apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->X, thSpecific->Y, thSpecific->costCache);
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                n2OptMoves = apply2OptBestFix_fastIteratively(&thSpecific->workingSol, thSpecific->costCache);
            #endif
        }

        if (nonImprovingIterCount > restartThreshold)
        {
            setupThSpecificOnBestSol(thSpecific);
            nonImprovingIterCount = 0;
        }
        if (thSpecific->workingSol.cost < thShared->bestSol->cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol->cost)
            {
                LOG(LOG_LVL_INFO, "Found better solution: cost = %lf", cvtCost2Double(thSpecific->workingSol.cost));
                cloneSolution(&thSpecific->workingSol, thShared->bestSol);
            }
            else
                nonImprovingIterCount++;
            pthread_mutex_unlock(&thShared->mutex);
        }
        else
            nonImprovingIterCount++;

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

    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
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

    for (int i = 0; i < n; i++) // build cost cache
    {
        #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            thSpecific->costCache[i] = computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[i+1], thSpecific->Y[i+1], inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            int *path = thShared->bestSol->indexPath;
            thSpecific->costCache[i] = inst->edgeCostMat[path[i] * n + path[i+1]];
        #endif
    }

    for (int i = n+1; i < n + AVX_VEC_SIZE; i++)
        thSpecific->costCache[i] = INFINITY;
    
    // reset tenure
    for (int i = 0; i < thShared->tenureSize; i++)
        thSpecific->tenure[i].node0 = -1;
    
    thSpecific->nextTenurePos = 0;
    thSpecific->lastTenurePos = 0;
}

static inline int randomlySelectEdgeOutsideTenure(ThreadSpecificData *thSpecific, int forbiddenNode)
{
    int edge = forbiddenNode, n = thSpecific->workingSol.instance->nNodes;

    // cannot choose any edge in the tenure and edge with the same value as 
    while ((edge == forbiddenNode) || (thSpecific->costCache[edge] == -INFINITY))
        edge = genRandom(&thSpecific->rndState, 0, n);

    return edge;
}

static void addEdgeToTenure(ThreadSpecificData *thSpecific, int indexPathIndex)
{
    #ifdef DEBUG
        if ((indexPathIndex < 0 ) || (indexPathIndex >= thSpecific->workingSol.instance->nNodes))
            throwError("Tabu -> addEdgeToTenure: indexPathIndex is not valid : %d", indexPathIndex);
    #endif

    if (thSpecific->nextTenurePos == thSpecific->lastTenurePos) // tenure full: remove last element from tenure and then add the new one
    {
        for (int i = 0; i < thSpecific->workingSol.instance->nNodes; i++)
        {
            if (thSpecific->workingSol.indexPath[i] == thSpecific->tenure[thSpecific->nextTenurePos].node0)
            {
                if ((thSpecific->tenure[thSpecific->nextTenurePos].node1 == thSpecific->workingSol.indexPath[i-1]) && (i > 0))
                    i--;
                thSpecific->costCache[i] = thSpecific->costBackup[thSpecific->nextTenurePos];
                break;
            }
        }
        thSpecific->lastTenurePos++;
        if (thSpecific->lastTenurePos == thSpecific->thShared->tenureSize)
            thSpecific->lastTenurePos = 0;
    }

    // add one element to the tenure
    thSpecific->tenure[thSpecific->nextTenurePos].node0 = thSpecific->workingSol.indexPath[indexPathIndex];
    thSpecific->tenure[thSpecific->nextTenurePos].node1 = thSpecific->workingSol.indexPath[indexPathIndex+1];
    thSpecific->costBackup[thSpecific->nextTenurePos] = thSpecific->costCache[indexPathIndex];
    thSpecific->costCache[indexPathIndex] = -INFINITY;

    thSpecific->nextTenurePos++;
    if (thSpecific->nextTenurePos == thSpecific->thShared->tenureSize)
        thSpecific->nextTenurePos = 0;
}

static void removeOldestElementFromTenure(ThreadSpecificData *thSpecific)
{   
    int n = thSpecific->workingSol.instance->nNodes;

    if (thSpecific->lastTenurePos == thSpecific->nextTenurePos)
        return; // tenure empty
    
    // replace correct cost into costCache
    for (int i = 0; i < n; i++)
    {
        if (thSpecific->workingSol.indexPath[i] == thSpecific->tenure[thSpecific->lastTenurePos].node0)
        {
            if ((thSpecific->tenure[thSpecific->lastTenurePos].node1 == thSpecific->workingSol.indexPath[i-1]) && (i > 0))
                i--;
            thSpecific->costCache[i] = thSpecific->costBackup[thSpecific->lastTenurePos];
            break;
        }
    }
    
    thSpecific->tenure[thSpecific->lastTenurePos].node0 = -1;

    thSpecific->lastTenurePos++;
    if (thSpecific->lastTenurePos == thSpecific->thShared->tenureSize)
            thSpecific->lastTenurePos = 0;
}

static inline void performNonImproving2OptMove(ThreadSpecificData *thSpecific, int edge0, int edge1)
{
    Solution *sol = &thSpecific->workingSol;
    Instance *inst = sol->instance;

    if (edge0 > edge1)
        swapElems(edge0, edge1)

    float altEdge0Cost, altEdge1Cost;

    #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        altEdge0Cost = computeEdgeCost(thSpecific->X[edge0], thSpecific->Y[edge0], thSpecific->X[edge1], thSpecific->Y[edge1], inst);
        altEdge1Cost = computeEdgeCost(thSpecific->X[edge0+1], thSpecific->Y[edge0+1], thSpecific->X[edge1+1], thSpecific->Y[edge1+1], inst);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        int *indexPath = sol->indexPath;
        altEdge0Cost = inst->edgeCostMat[(size_t)indexPath[edge0] * (size_t)inst->nNodes + (size_t)indexPath[edge1]];
        altEdge1Cost = inst->edgeCostMat[(size_t)indexPath[edge0+1] * (size_t)inst->nNodes + (size_t)indexPath[edge1+1]];
    #endif

    LOG(LOG_LVL_TRACE, "Tabu[%d]: Updating solution by switching edge (%d,%d) with edge (%d,%d) degrading cost by %f. New Cost = %lf", thSpecific->iterCount,
        thSpecific->workingSol.indexPath[edge0], thSpecific->workingSol.indexPath[edge0+1],
        thSpecific->workingSol.indexPath[edge1], thSpecific->workingSol.indexPath[edge1+1],
        altEdge0Cost + altEdge1Cost - thSpecific->costCache[edge0] - thSpecific->costCache[edge1], cvtCost2Double(sol->cost));

    // update cost
    sol->cost += cvtFloat2Cost(altEdge0Cost) + cvtFloat2Cost(altEdge1Cost) - cvtFloat2Cost(thSpecific->costCache[edge0]) - cvtFloat2Cost(thSpecific->costCache[edge1]);

    int smallID = edge0 + 1, bigID = edge1;

    while (smallID < bigID)
    {
        swapElems(sol->indexPath[smallID], sol->indexPath[bigID])
        #if ((COMPUTATION_TYPE == COMPUTATE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            swapElems(thSpecific->X[smallID], thSpecific->X[bigID])
            swapElems(thSpecific->Y[smallID], thSpecific->Y[bigID])
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
        swapElems(thSpecific->costCache[smallID], thSpecific->costCache[bigID])

        smallID++;
        bigID--;
    }
}

#ifdef DEBUG
static void checkThSpecificData(ThreadSpecificData *thSpecific)
{
    Solution *sol = &thSpecific->workingSol;
    Instance *inst = sol->instance;

    // check tenure
    for (int i = thSpecific->lastTenurePos; i != thSpecific->nextTenurePos; i++)
    {
        if (i >= thSpecific->thShared->tenureSize)
        {
            i = -1;
            continue;
        }
        if ((thSpecific->tenure[i].node0 < 0) || (thSpecific->tenure[i].node0 >= inst->nNodes))
            throwError("Tabu -> Tenure has node0 ad position %d set to %d which is invalid", i, thSpecific->tenure[i].node0);
        float recompCost = computeEdgeCost(inst->X[thSpecific->tenure[i].node0], inst->Y[thSpecific->tenure[i].node0], inst->X[thSpecific->tenure[i].node1], inst->Y[thSpecific->tenure[i].node1], inst);
        if (thSpecific->costBackup[i] != recompCost)
            throwError("Tabu -> CostBackup contains %f while it should contain %f", thSpecific->costBackup[i], recompCost);
    }
    
    for (int i = 0; i < inst->nNodes; i++)
    {
        bool notInTenure = true;
        for (int j = 0; j < thSpecific->thShared->tenureSize; j++)
        {
            if (thSpecific->tenure[j].node0 == -1) continue;

            if (((sol->indexPath[i] == thSpecific->tenure[j].node0) && (sol->indexPath[i+1] == thSpecific->tenure[j].node1)) ||
                ((sol->indexPath[i] == thSpecific->tenure[j].node1) && (sol->indexPath[i+1] == thSpecific->tenure[j].node0)))
            {
                if (thSpecific->costCache[i] != -INFINITY)
                    throwError("Tabu -> checkTenureAndLocks: Edge[%d,%d] in tenure is NOT set to -INFINITY inside the costCache", sol->indexPath[i], sol->indexPath[i+1]);
                notInTenure = false;
            }
        }
        if (notInTenure)
        {
            if (thSpecific->costCache[i] == -INFINITY)
                throwError("Tabu -> checkTenureAndLocks: Edge[%d,%d] NOT in tenure is set to -INFINITY inside the costCache", sol->indexPath[i], sol->indexPath[i+1]);
            else
            {
                float recomputedCost;
                #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                    recomputedCost = computeEdgeCost(thSpecific->X[i], thSpecific->Y[i], thSpecific->X[i+1], thSpecific->Y[i+1], inst);
                #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                    recomputedCost = inst->edgeCostMat[sol->indexPath[i] * inst->nNodes + sol->indexPath[i+1]];
                #endif
                if (thSpecific->costCache[i] != recomputedCost)
                    throwError("Tabu -> checkTenureAndLocks: Cost cache is not coherent costCache[%d] = %f which is not %f", i, thSpecific->costCache[i], recomputedCost);
            }
        }
    }
}

static void printTenure(ThreadSpecificData *thSpecific, char *strOut)
{
    if (thSpecific->workingSol.instance->params.logLevel < LOG_LVL_DEBUG)
        return;

    int tenureSize = thSpecific->thShared->tenureSize;
    sprintf(strOut, "Tenure: [");
    strOut = &strOut[9];
    bool tenureIsEmpty = true;

    for (int i = thSpecific->nextTenurePos; i < tenureSize; i++)
        if (thSpecific->tenure[i].node0 != -1)
        {
            sprintf(strOut, "(%2d,%2d), ", thSpecific->tenure[i].node0, thSpecific->tenure[i].node1);
            strOut = &strOut[9];
            tenureIsEmpty = false;
        }

    for (int i = 0; i < thSpecific->nextTenurePos; i++)
        if (thSpecific->tenure[i].node0 != -1)
        {
            sprintf(strOut, "(%2d,%2d), ", thSpecific->tenure[i].node0, thSpecific->tenure[i].node1);
            strOut = &strOut[9];
            tenureIsEmpty = false;
        }

    if (tenureIsEmpty)
        sprintf(strOut, "]");
    else
        sprintf(strOut, "\b\b]");
}
#endif
