#include "Tsp.h"


#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

typedef struct
{
    Solution bestSol;

    enum EMInitType startOption;

    pthread_mutex_t mutex;

    double timeLimit;
} ThreadSharedData;

static ThreadSharedData initThreadSharedData (Instance *inst, enum EMInitType startOption, double timeLimit)
{
    ThreadSharedData thShared = {
        .startOption = startOption,
        .timeLimit = timeLimit
    };

    if (pthread_mutex_init(&thShared.mutex, NULL)) throwError(inst, NULL, "ExtraMileage -> initThreadSharedData: Failed to initialize mutex");

    thShared.bestSol = newSolution(inst);

    return thShared;
}

static void destroyThreadSharedData (ThreadSharedData *thShared)
{
    if (pthread_mutex_init(&thShared->mutex, NULL)) throwError(thShared->bestSol.instance, &thShared->bestSol, "ExtraMileage -> destroyThreadSharedData: Failed to destroy mutex");
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

    thSpecific.X = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(int));
    if (!thSpecific.X)
        throwError(inst, &thSpecific.workingSol, "ExtraMileage -> initThreadSpecificData: Failed to allocate memory");
    thSpecific.Y = &thSpecific.X[inst->nNodes + AVX_VEC_SIZE];

    return thSpecific;
}

static void destroyThreadSpecificData(ThreadSpecificData *thSpecific)
{
    destroySolution(&thSpecific->workingSol);

    free(thSpecific->X);
    thSpecific->X = thSpecific->Y = NULL;
}

static inline void swapElemsInThSpecific(ThreadSpecificData *thSpecific, int pos1, int pos2)
{
    register float tempFloat;
    swapElems(thSpecific->X[pos1], thSpecific->X[pos2], tempFloat);
    swapElems(thSpecific->Y[pos1], thSpecific->Y[pos2], tempFloat);
    register int tempInt;
    swapElems(thSpecific->workingSol.indexPath[pos1], thSpecific->workingSol.indexPath[pos2], tempInt);
}


static void *runExtraMileage(void *arg);

Solution ExtraMileage(Instance *inst, enum EMInitType startOption, double timeLimit, int nThreads)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    /*if ((COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX) && (inst->edgeCostMat == NULL))
    {
        LOG(LOG_LVL_WARNING, "ExtraMileage: edgeCostMat has not been detected. Computing it now");
        double costMatTime = computeCostMatrix(inst);
        LOG(LOG_LVL_WARNING, "ExtraMileage: edgeCostMat computation took: %lf seconds", costMatTime);
    }*/

    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError(inst, NULL, "ExtraMileage: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = inst->params.nThreads;

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


static void applyExtraMileage_Internal(ThreadSpecificData *thSpecific, int nCovered);

void applyExtraMileage(Solution *sol, int nCovered, unsigned int *rndState)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    // check input solution
    for (int i = 0; i < nCovered; i++)
        for (int j = 0; j < nCovered; j++)
            if ((i != j) && (sol->indexPath[i] == sol->indexPath[j]))
                throwError(sol->instance, sol, "applyExtraMileage: Presented solution has the same node in the path before nCovered. sol->indexPath[%d] = sol->indexPath[%d] = %d", i, j, sol->indexPath[i]);

    // recompute cost of nCovered loop
    sol->cost = 0;
    for (int i = 0; i < nCovered-1; i++)
        sol->cost += computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i+1]], inst->Y[sol->indexPath[i+1]], inst->params.edgeWeightType, inst->params.roundWeights);
    sol->cost += computeEdgeCost(inst->X[sol->indexPath[nCovered-1]], inst->Y[sol->indexPath[nCovered-1]], inst->X[sol->indexPath[0]], inst->Y[sol->indexPath[0]], inst->params.edgeWeightType, inst->params.roundWeights);
            
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

    cloneSolution(sol, &thSpecific.workingSol);
    for (int i = 0; i < (n + AVX_VEC_SIZE) * 2; i++)
        thSpecific.X[i] = inst->X[sol->indexPath[i]];

    applyExtraMileage_Internal(&thSpecific, nCovered);

    cloneSolution(&thSpecific.workingSol, sol);

    *rndState = thSpecific.rndState;
    destroyThreadSpecificData(&thSpecific);
}


static int initialization(ThreadSpecificData *thSpecific);

static void updateBestSolution(ThreadSpecificData *thSpecific);

static void *runExtraMileage(void * arg)
{
    ThreadSpecificData *thSpecific = (ThreadSpecificData*)arg;
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;

    // optimization: save initialization if it is the same for every run (if init mode is random then it must be generated every time so no optimization)
    Solution init;
    int nCovered;
    if (thShared->startOption != EM_INIT_RANDOM)
    {
        init = newSolution(inst);
        nCovered = initialization(thSpecific);
    }

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    while (currentTime < thShared->timeLimit)
    {
        if (thShared->startOption != EM_INIT_RANDOM)
            cloneSolution(&init, &thSpecific->workingSol);
        else
            nCovered = initialization(thSpecific);
        
        applyExtraMileage_Internal(thSpecific, nCovered);

        if (thSpecific->workingSol.cost < thShared->bestSol.cost)
        {
            pthread_mutex_lock(&thShared->mutex);
            if (thSpecific->workingSol.cost < thShared->bestSol.cost)
                updateBestSolution(thSpecific);
            pthread_mutex_unlock(&thShared->mutex);
        }

        thSpecific->iterCount++;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }

    return NULL;
}

static void farthestPointsInit(ThreadSpecificData *thSpecific);

static int initialization(ThreadSpecificData *thSpecific)
{
    ThreadSharedData *thShared = thSpecific->thShared;
    Instance *inst = thSpecific->workingSol.instance;
    
    //thSpecific->workingSol.cost = 0;

    // Set data for solution (copy coords from distance and create index path as 0,1,2,...,n-1)
    for (int i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
        thSpecific->X[i] = inst->X[i];
    for (int i = 0; i < inst->nNodes + AVX_VEC_SIZE; i++)
        thSpecific->workingSol.indexPath[i] = i;

    int nCovered;

    switch (thShared->startOption)
    {
    case EM_INIT_RANDOM:
        {
        // select two random nodes
        int rndIndex0 = genRandom(&thSpecific->rndState, 0, inst->nNodes), rndIndex1 = genRandom(&thSpecific->rndState, 0, inst->nNodes);
        while (abs(rndIndex1 - rndIndex0) <= 1)
            rndIndex1 = genRandom(&thSpecific->rndState, 0, inst->nNodes);

        // add the two nodes to the solution (order does not matter)
        swapElemsInThSpecific(thSpecific, 0, rndIndex0);
        swapElemsInThSpecific(thSpecific, 1, rndIndex1);

        // update cost
        thSpecific->workingSol.cost = computeEdgeCost(thSpecific->X[0], thSpecific->Y[0], thSpecific->X[1], thSpecific->Y[1], inst->params.edgeWeightType , inst->params.roundWeights) * 2.;

        //LOG(LOG_LVL_DEBUG, "ExtraMileage-InitializationRandom: Randomly chosen edge is edge e = (%d, %d)", rndIndex0, rndIndex1);

        nCovered = 2;
        }
        break;

    case EM_INIT_FARTHEST_POINTS:
        farthestPointsInit(thSpecific);

        // update cost
        thSpecific->workingSol.cost = computeEdgeCost(thSpecific->X[0], thSpecific->Y[0], thSpecific->X[1], thSpecific->Y[1], inst->params.edgeWeightType , inst->params.roundWeights) * 2.;

        nCovered = 2;
        break;
    }

    return nCovered;
}

static void farthestPointsInit(ThreadSpecificData *thSpecific)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    __m256 maxCostVec = _mm256_set1_ps(0), rowMaxCostVec = _mm256_set1_ps(0); // cost is always positive
    __m256i maxIndexVec1 = _mm256_set1_epi32(0), maxIndexVec2 = _mm256_set1_epi32(0);
    __m256i incrementVec = _mm256_set1_epi32(AVX_VEC_SIZE), ones = _mm256_set1_epi32(1);
    

    __m256i rowIDsVec = _mm256_set1_epi32(0); // the content of this are always all i
    for (int i = 0; i < n-1; i++, rowIDsVec = _mm256_add_epi32(rowIDsVec, ones))
    {
        __m256i colIDsVec = _mm256_set_epi32( 8+i, 7+i, 6+i, 5+i, 4+i, 3+i, 2+i, 1+i );
        __m256 x1 = _mm256_broadcast_ss(&inst->X[i]), y1 = _mm256_broadcast_ss(&inst->Y[i]);
        
        for (int j = i+1; j < n; j += AVX_VEC_SIZE, colIDsVec = _mm256_add_epi32(colIDsVec, incrementVec))
        {
            if (j > n - AVX_VEC_SIZE)
            {
                colIDsVec = _mm256_sub_epi32(colIDsVec, _mm256_set1_epi32(AVX_VEC_SIZE - n % AVX_VEC_SIZE));
                j = n - AVX_VEC_SIZE; // subtract AVX_VEC_SIZE one extra time to compensate the increment of the loop
            }

            __m256 x2 = _mm256_loadu_ps(&inst->X[j]), y2 = _mm256_loadu_ps(&inst->Y[j]);
            __m256 costVec = computeEdgeCost_VEC(x1, y1, x2, y2, inst->params.edgeWeightType , inst->params.roundWeights);

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
    _mm256_storeu_si256((__m256i_u*)maxIndexes, maxIndexVec1);
    int maxIndex1 = maxIndexes[maxIndex];
    _mm256_storeu_si256((__m256i_u*)maxIndexes, maxIndexVec2);
    int maxIndex2 = maxIndexes[maxIndex];

    LOG(LOG_LVL_DEBUG, "Extra Mileage EM_INIT_FARTHEST_POINTS: maximum cost found is %f between nodes %d and %d", *(float*)&maxCosts[maxIndex], maxIndex1, maxIndex2);

    // initialize solution
    swapElemsInThSpecific(thSpecific, 0, maxIndex1);
    swapElemsInThSpecific(thSpecific, 1, maxIndex2);
}

static void updateBestSolution(ThreadSpecificData *thSpecific)
{
    Solution *bestSol = &thSpecific->thShared->bestSol;
    Solution *newBest = &thSpecific->workingSol;

    // check solution when debugging
    if (bestSol->instance->params.logLevel >= LOG_LVL_DEBUG)
        if (!checkSolution(newBest))
        {
            destroySolution(bestSol);
		    throwError(newBest->instance, newBest, "updateBestSolutionEM: newBest Solution is not valid");
        }
    
    LOG(LOG_LVL_LOG, "Found better solution: cost = %f", newBest->cost);

    bestSol->cost = newBest->cost;

    register int *temp;
    swapElems(bestSol->indexPath, newBest->indexPath, temp);
}


static bool checkSolutionIntegrity(ThreadSpecificData *thSpecific);

typedef struct 
{
    int node;
    int anchor;
    float extraCost;
} SuccessorData;

static inline void insertNodeInSolution(ThreadSpecificData *thSpecific, int nCovered, SuccessorData succ);

static SuccessorData findSuccessorVectorized(ThreadSpecificData *thSpecific, int nCovered);

static SuccessorData findSuccessorBase(ThreadSpecificData *thSpecific, int nCovered, bool useCostMatrix);

static void applyExtraMileage_Internal(ThreadSpecificData *thSpecific, int nCovered)
{
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;

    if (!checkSolutionIntegrity(thSpecific))
        throwError(inst, &thSpecific->workingSol, "applyExtraMileage_Internal: thSpecific.workingSol is not consistent");

    // save element to last position
    thSpecific->X[n] = thSpecific->X[nCovered];
    thSpecific->Y[n] = thSpecific->Y[nCovered];
    thSpecific->workingSol.indexPath[n] = thSpecific->workingSol.indexPath[nCovered];

    // close the tour at index nCovered
    thSpecific->X[nCovered] = thSpecific->X[0];
    thSpecific->Y[nCovered] = thSpecific->Y[0];
    thSpecific->workingSol.indexPath[nCovered] = thSpecific->workingSol.indexPath[0];

    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    for (; nCovered < n; nCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        SuccessorData succ;

        if ((inst->params.graspType == GRASP_RANDOM) && (graspThreshold > rand_r(&thSpecific->rndState)))
        {
            enum EdgeWeightType ewt = thSpecific->workingSol.instance->params.edgeWeightType;
            bool roundW = thSpecific->workingSol.instance->params.roundWeights;
            succ.node = genRandom(&thSpecific->rndState, (nCovered + 1), (n + 1));
            succ.anchor = genRandom(&thSpecific->rndState, 0, nCovered);
            succ.extraCost = computeEdgeCost(thSpecific->X[succ.anchor], thSpecific->Y[succ.anchor], thSpecific->X[succ.node],       thSpecific->Y[succ.node],       ewt, roundW) +
                             computeEdgeCost(thSpecific->X[succ.node],   thSpecific->Y[succ.node],   thSpecific->X[succ.anchor + 1], thSpecific->Y[succ.anchor + 1], ewt, roundW) -
                             computeEdgeCost(thSpecific->X[succ.anchor], thSpecific->Y[succ.anchor], thSpecific->X[succ.anchor + 1], thSpecific->Y[succ.anchor + 1], ewt, roundW);
        }
        else
        {
            switch (COMPUTATION_TYPE)
            {
            case COMPUTE_OPTION_AVX:
                succ = findSuccessorVectorized(thSpecific, nCovered);
                break;
            case COMPUTE_OPTION_BASE:
                succ = findSuccessorBase(thSpecific, nCovered, false);
                break;
            case COMPUTE_OPTION_USE_COST_MATRIX:
                succ = findSuccessorBase(thSpecific, nCovered, true);
                break;
            }
        }

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
        if ((thSpecific->X[i] != inst->X[index]) || (thSpecific->Y[i] != inst->Y[index]))
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: Mismatch at index %lu in solution", i);
            return false;
        }
    }

    // everything checks out
    return true;
}

static inline void insertNodeInSolution(ThreadSpecificData *thSpecific, int nCovered, SuccessorData succ)
{
    nCovered++;

    // update cost
    thSpecific->workingSol.cost += succ.extraCost;

    // save best value
    float bestX = thSpecific->X[succ.node];
    float bestY = thSpecific->Y[succ.node];
    int bestIndex = thSpecific->workingSol.indexPath[succ.node];

    // place elements to insert in the tour at the end of the covered nodes "set"
    thSpecific->X[succ.node] = thSpecific->X[nCovered];
    thSpecific->Y[succ.node] = thSpecific->Y[nCovered];
    thSpecific->workingSol.indexPath[succ.node] = thSpecific->workingSol.indexPath[nCovered];

    int i = nCovered;

    // shift elements forward of 1 position iteratively with avx until vector is too big for the amount of elements to shift (do AVX_VEC_SIZE elements per iteration)
    for (i -= AVX_VEC_SIZE; i > succ.anchor; i -= AVX_VEC_SIZE)
    {
        __m256 xVec = _mm256_loadu_ps(&thSpecific->X[i]);
        __m256 yVec = _mm256_loadu_ps(&thSpecific->Y[i]);
        __m256i indexVec = _mm256_loadu_si256((__m256i_u*)&thSpecific->workingSol.indexPath[i]);
        _mm256_storeu_ps(&thSpecific->X[i + 1], xVec);
        _mm256_storeu_ps(&thSpecific->Y[i + 1], yVec);
        _mm256_storeu_si256((__m256i_u*)&thSpecific->workingSol.indexPath[i + 1], indexVec);
    }
    i += AVX_VEC_SIZE;

    // shift elements forward one at a time
    for (i--; i > succ.anchor; i--)
    {
        thSpecific->X[i+1] = thSpecific->X[i];
        thSpecific->Y[i+1] = thSpecific->Y[i];
        thSpecific->workingSol.indexPath[i+1] = thSpecific->workingSol.indexPath[i];
    }

    i++;

    thSpecific->X[i] = bestX;
    thSpecific->Y[i] = bestY;
    thSpecific->workingSol.indexPath[i] = bestIndex;

    LOG(LOG_LVL_EVERYTHING, "Extra Mileage Solution Update: Node %d added to solution between nodes %d and %d", thSpecific->workingSol.indexPath[i], thSpecific->workingSol.indexPath[i-1], thSpecific->workingSol.indexPath[i+1]);
}

static SuccessorData findSuccessorVectorized(ThreadSpecificData *thSpecific, int nCovered)
{
    // shortcuts/decluttering
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundW = inst->params.roundWeights;
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    // Contains best mileage values
    __m256 bestExtraMileageVec = _mm256_set1_ps(INFINITY);
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
        __m256 curEdgeCostVec = computeEdgeCost_VEC(x1Vec, y1Vec, x2Vec, y2Vec, ewt, roundW);

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
            __m256 cmpMask = _mm256_cmp_ps(curExtraMileageVec, bestExtraMileageVec, _CMP_LT_OQ);

            // Set new best according to comparison result
            bestExtraMileageVec = _mm256_blendv_ps(bestExtraMileageVec, curExtraMileageVec, cmpMask);
            bestAnchorsVec = _mm256_blendv_epi8(bestAnchorsVec, curEdgeID, _mm256_castps_si256(cmpMask));
            bestNodesVec = _mm256_blendv_epi8(bestNodesVec, idsVec, _mm256_castps_si256(cmpMask));
        }
    }
    // at this point we must select a candidate(best or almost-best(grasp))

    float bestExtraMileage[AVX_VEC_SIZE];
    int bestNodes[AVX_VEC_SIZE];
    int bestAnchors[AVX_VEC_SIZE];

    _mm256_storeu_ps(bestExtraMileage, bestExtraMileageVec);
    _mm256_storeu_si256((__m256i_u *)bestNodes, bestNodesVec);
    _mm256_storeu_si256((__m256i_u *)bestAnchors, bestAnchorsVec);

    int sortedArgs[AVX_VEC_SIZE];
    argsort(bestExtraMileage, sortedArgs, AVX_VEC_SIZE);

    int chosedNodeSubIndex = 0;
    if ((inst->params.graspType == GRASP_ALMOSTBEST) && (graspThreshold > rand_r(&thSpecific->rndState)) && (n - nCovered - 1 > AVX_VEC_SIZE))
    {
        chosedNodeSubIndex = 1;
        for (; chosedNodeSubIndex < AVX_VEC_SIZE - 1; chosedNodeSubIndex++)
            if (rand_r(&thSpecific->rndState) > RAND_MAX / 2)
                break;
    }

    SuccessorData retVal = { .node = bestNodes[sortedArgs[chosedNodeSubIndex]], .anchor = bestAnchors[sortedArgs[chosedNodeSubIndex]], .extraCost = bestExtraMileage[sortedArgs[chosedNodeSubIndex]] };

    return retVal;
}

#define BASE_GRASP_BEST_SAVE_BUFFER_SIZE 4
static SuccessorData findSuccessorBase(ThreadSpecificData *thSpecific, int nCovered, bool useCostMatrix)
{
    // shortcuts/decluttering
    Instance *inst = thSpecific->workingSol.instance;
    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundW = inst->params.roundWeights;
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    SuccessorData bestSuccs[BASE_GRASP_BEST_SAVE_BUFFER_SIZE];
    for (int i = 0; i < BASE_GRASP_BEST_SAVE_BUFFER_SIZE; i++)
        bestSuccs[i].extraCost = INFINITY;
    
    for (int anchor = 0; anchor < nCovered; anchor++) // u stands for uncovered
    {
        // cost of edge already in solution [i,j]
        float currEdgeCost;
        if (useCostMatrix)
            currEdgeCost = inst->edgeCostMat[thSpecific->workingSol.indexPath[anchor] * n + thSpecific->workingSol.indexPath[anchor + 1]];
        else
            currEdgeCost = computeEdgeCost(thSpecific->X[anchor], thSpecific->Y[anchor], thSpecific->X[anchor + 1], thSpecific->Y[anchor + 1], ewt, roundW);

        for (int node = nCovered+1; node < n + 1; node++) // covered node i from 0 to posCovered
        {
            // sum of cost of edge [i,u] and edge [u,j]
            float altEdgeCost;
            if (useCostMatrix)
                altEdgeCost = inst->edgeCostMat[thSpecific->workingSol.indexPath[node] * n + thSpecific->workingSol.indexPath[anchor]] + 
                              inst->edgeCostMat[thSpecific->workingSol.indexPath[node] * n + thSpecific->workingSol.indexPath[anchor + 1]];
            else
                altEdgeCost = computeEdgeCost(thSpecific->X[node], thSpecific->Y[node], thSpecific->X[anchor], thSpecific->Y[anchor], ewt, roundW) + 
                              computeEdgeCost(thSpecific->X[node], thSpecific->Y[node], thSpecific->X[anchor + 1], thSpecific->Y[anchor + 1], ewt, roundW);

            SuccessorData currSucc = { .node=node, .anchor=anchor, .extraCost = altEdgeCost - currEdgeCost };

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
