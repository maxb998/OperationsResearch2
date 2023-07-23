#include "Tsp.h"


#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// Flag to define at compile time whether to use or not the AVX instruction to move the data
#define EM_USE_FAST_SOLUTION_UPDATE 1

typedef struct
{
    Solution *bestSol;

    enum EMInitType emInitType;
    enum EMOptions emOpt;

    pthread_mutex_t mutex;

    double tlim;
} EMThreadsSharedData;

typedef struct
{
    EMThreadsSharedData *shared;

    unsigned int *rndState;
    unsigned long iterCount;
} EMThreadsData;

static void *runExtraMileage(void *arg);

static void updateBestSolutionEM(Solution *bestSol, Solution *newBest);

static size_t initialization(Solution *sol, enum EMInitType emType, unsigned int *rndState);

static void farthestPointsInit(Solution *sol);

static bool checkSolutionIntegrity(Solution *sol);

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2);

static inline void updateSolutionEM(Solution *sol, size_t posCovered, size_t bestMileageIndex, size_t anchorIndex, float extraCost);

static void extraMileageVectorized(Solution *sol, size_t nCovered, unsigned int *rndState);

static void extraMileageBase(Solution *sol, size_t nCovered, unsigned int *rndState, bool useCostMatrix);


Solution ExtraMileage(Instance *inst, enum EMOptions emOpt, enum EMInitType emInitType, double tlim, int nThreads)
{
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError(inst, NULL, "ExtraMileage: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = inst->params.nThreads;

    // apply extra mileage
    Solution bestSol = newSolution(inst);
    EMThreadsSharedData shared = { .bestSol=&bestSol, .emInitType=emInitType, .emOpt=emOpt, .tlim=tlim };
    unsigned long iterCount = 0;

    pthread_mutex_init(&shared.mutex, NULL);
    EMThreadsData th[MAX_THREADS];
    unsigned int states[MAX_THREADS];
    for (size_t i = 0; i < nThreads; i++)
    {
        th[i].shared = &shared;
        states[i] = (unsigned int)rand();
        th[i].rndState = &states[i];
        th[i].iterCount = 0;
    }

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < nThreads; i++)
        pthread_create(&threads[i], NULL, runExtraMileage, &th[i]);

    for (size_t i = 0; i < nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        iterCount += th[i].iterCount;
    }

    pthread_mutex_destroy(&shared.mutex);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    bestSol.execTime = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);

    LOG(LOG_LVL_NOTICE, "Total number of iterations: %ld", iterCount);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterCount/bestSol.execTime);

    return bestSol;
}

static void *runExtraMileage(void * arg)
{
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double tStart = currT.tv_sec + currT.tv_nsec / 1000000000.0;
    double t = tStart;

    EMThreadsData *th = (EMThreadsData*)arg;
    EMThreadsSharedData *shared = th->shared;
    Solution *bestSol = shared->bestSol;
    Instance *inst = bestSol->instance;

    // optimization: save initialization if it is the same for every run (if init mode is random then it must be generated every time so no optimization)
    Solution init;
    size_t coveredNodes;
    if (shared->emInitType != EM_INIT_RANDOM)
    {
        init = newSolution(inst);
        coveredNodes = initialization(&init, shared->emInitType, th->rndState);
    }

    Solution sol = newSolution(inst);

    while (t < tStart + shared->tlim)
    {
        if (inst->params.emInitOption != EM_INIT_RANDOM)
            cloneSolution(&init, &sol);
        else
            coveredNodes = initialization(&sol, EM_INIT_RANDOM, th->rndState);
        
        applyExtraMileage(&sol, coveredNodes, shared->emOpt, th->rndState);

        if (sol.cost < bestSol->cost)
        {
            pthread_mutex_lock(&shared->mutex);
            if (sol.cost < bestSol->cost)
                updateBestSolutionEM(bestSol, &sol);
            pthread_mutex_unlock(&shared->mutex);
        }

        th->iterCount++;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        t = currT.tv_sec + currT.tv_nsec / 1000000000.0;
    }

    pthread_exit(NULL);
}

static void updateBestSolutionEM(Solution *bestSol, Solution *newBest)
{
    // check solution when debugging
    if (bestSol->instance->params.logLevel >= LOG_LVL_EVERYTHING)
        if (!checkSolution(newBest))
        {
            destroySolution(bestSol);
		    throwError(newBest->instance, newBest, "updateBestSolutionEM: newBest Solution is not valid");
        }
    
    LOG(LOG_LVL_LOG, "Found better solution: New cost: %f   Old cost: %f", newBest->cost, bestSol->cost);

    bestSol->cost = newBest->cost;

    register float *tempf;
    swapElems(bestSol->X, newBest->X, tempf);
    swapElems(bestSol->Y, newBest->Y, tempf);
    register int *tempi;
    swapElems(bestSol->indexPath, newBest->indexPath, tempi);
}

void applyExtraMileage(Solution *sol, size_t nCovered, enum EMOptions emOpt, unsigned int *rndState)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    // first check that the solution element at index nCovered is the same as the one at index 0
    if (sol->indexPath[0] != sol->indexPath[nCovered])
    {
        // check solution integrity when debugging
        if (inst->params.logLevel >= LOG_LVL_DEBUG || !(sol->X[n] == INFINITY && sol->Y[n] == INFINITY) || !(sol->X[n] == sol->X[0] && sol->Y[n] == sol->Y[0]))
            if (!checkSolutionIntegrity(sol))
                throwError(inst, sol, "applyExtraMileage: Error when checking input solution");

        // save element to last position
        sol->X[n] = sol->X[nCovered];
        sol->Y[n] = sol->Y[nCovered];
        sol->indexPath[n] = sol->indexPath[nCovered];

        // close the tour at index nCovered
        sol->X[nCovered] = sol->X[0];
        sol->Y[nCovered] = sol->Y[0];
        sol->indexPath[nCovered] = sol->indexPath[0];
    }
    
    switch (emOpt)
    {
    case EM_OPTION_AVX:
        extraMileageVectorized(sol, nCovered, rndState);
        break;
    case EM_OPTION_BASE:
        extraMileageBase(sol, nCovered, rndState, false);
        break;
    case EM_OPTION_USE_COST_MATRIX:
        if (inst->edgeCostMat)
            extraMileageBase(sol, nCovered, rndState, true);
        else
        {
            LOG(LOG_LVL_WARNING, "applyExtraMileage: edgeCostMat has not been detected/initialized. Switching to option: EM_OPTION_AVX");
            extraMileageVectorized(sol, nCovered, rndState);
        }
        break;
    }
}

static size_t initialization(Solution *sol, enum EMInitType emInitType, unsigned int *rndState)
{
    Instance *inst = sol->instance;
    size_t coveredElems = 0;
    sol->cost = 0;

    // Set data for solution (copy coords from distance and create index path as 0,1,2,...,n-1)
    for (size_t i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
        sol->X[i] = inst->X[i];
    for (int i = 0; i < inst->nNodes; i++)
        sol->indexPath[i] = i;

    switch (emInitType)
    {
    case EM_INIT_RANDOM:
        {
        // select two random nodes
        int rndIndex0 = rand_r(rndState) % (int)inst->nNodes, rndIndex1 = rand_r(rndState) % (int)inst->nNodes;
        while (abs(rndIndex1 - rndIndex0) <= 1)
            rndIndex1 = rand_r(rndState) % (int)inst->nNodes;

        // add the two nodes to the solution (order does not matter)
        swapElementsInSolution(sol, 0, rndIndex0);
        swapElementsInSolution(sol, 1, rndIndex1);

        // update cost
        sol->cost = computeEdgeCost(sol->X[0], sol->Y[0], sol->X[1], sol->Y[1], inst->params.edgeWeightType , inst->params.roundWeights) * 2.;

        //LOG(LOG_LVL_DEBUG, "ExtraMileage-InitializationRandom: Randomly chosen edge is edge e = (%d, %d)", rndIndex0, rndIndex1);

        coveredElems = 2;
        }
        break;

    case EM_INIT_FARTHEST_POINTS:
        farthestPointsInit(sol);

        // update cost
        sol->cost = computeEdgeCost(sol->X[0], sol->Y[0], sol->X[1], sol->Y[1], inst->params.edgeWeightType , inst->params.roundWeights) * 2.;

        coveredElems = 2;
        break;
    case EM_INIT_HULL:
        // Not yet supported
        break;
    }

    // close the tour
    sol->X[inst->nNodes] = sol->X[2];
    sol->Y[inst->nNodes] = sol->Y[2];
    sol->indexPath[inst->nNodes] = sol->indexPath[2];
    sol->X[2] = sol->X[0];
    sol->Y[2] = sol->Y[0];
    sol->indexPath[2] = sol->indexPath[0];

    return coveredElems;
}

static void farthestPointsInit(Solution *sol)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    __m256 maxCostVec = _mm256_set1_ps(0), rowMaxCostVec = _mm256_set1_ps(0); // cost is always positive
    __m256i maxIndexVec1 = _mm256_set1_epi32(0), maxIndexVec2 = _mm256_set1_epi32(0);
    __m256i incrementVec = _mm256_set1_epi32(AVX_VEC_SIZE), ones = _mm256_set1_epi32(1);
    

    __m256i rowIDsVec = _mm256_set1_epi32(0); // the content of this are always all i
    for (size_t i = 0; i < n-1; i++, rowIDsVec = _mm256_add_epi32(rowIDsVec, ones))
    {
        __m256i colIDsVec = _mm256_set_epi32( 8+i, 7+i, 6+i, 5+i, 4+i, 3+i, 2+i, 1+i );
        __m256 x1 = _mm256_broadcast_ss(&inst->X[i]), y1 = _mm256_broadcast_ss(&inst->Y[i]);
        
        for (size_t j = i+1; j < n; j += AVX_VEC_SIZE, colIDsVec = _mm256_add_epi32(colIDsVec, incrementVec))
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
    size_t maxIndex = 0;
    for (size_t i = 1; i < AVX_VEC_SIZE; i++)
        if (maxCosts[maxIndex] < maxCosts[i])
            maxIndex = i;

    int maxIndexes[AVX_VEC_SIZE];
    _mm256_storeu_si256((__m256i_u*)maxIndexes, maxIndexVec1);
    int maxIndex1 = maxIndexes[maxIndex];
    _mm256_storeu_si256((__m256i_u*)maxIndexes, maxIndexVec2);
    int maxIndex2 = maxIndexes[maxIndex];

    LOG(LOG_LVL_DEBUG, "Extra Mileage EM_INIT_FARTHEST_POINTS: maximum cost found is %f between nodes %d and %d", *(float*)&maxCosts[maxIndex], maxIndex1, maxIndex2);

    // initialize solution
    swapElementsInSolution(sol, 0, maxIndex1);
    swapElementsInSolution(sol, 1, maxIndex2);
}

static bool checkSolutionIntegrity(Solution *sol)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    for (size_t i = 0; i < n; i++)
    {
        int index = sol->indexPath[i];
        if (index < 0 || index > n)
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: sol.indexPath[%lu] = %d which is not within the limits", i, index);
            return false;
        }
        if (sol->X[i] != inst->X[index] || sol->Y[i] != inst->Y[index])
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: Mismatch at index %lu in solution", i);
            return false;
        }
    }

    // everything checks out
    return true;
}

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2)
{
    register float tempf;
    swapElems(sol->X[pos1], sol->X[pos2], tempf);
    swapElems(sol->Y[pos1], sol->Y[pos2], tempf);
    register int tempi;
    swapElems(sol->indexPath[pos1], sol->indexPath[pos2], tempi);
}

static inline void updateSolutionEM(Solution *sol, size_t posCovered, size_t bestMileageIndex, size_t anchorIndex, float extraCost)
{
    // update cost
    sol->cost += extraCost;

    // save best value
    float bestX = sol->X[bestMileageIndex], bestY = sol->Y[bestMileageIndex];
    int bestIndex = sol->indexPath[bestMileageIndex];

    // place elements to insert in the tour at the end of the covered nodes "set"
    sol->X[bestMileageIndex] = sol->X[posCovered];
    sol->Y[bestMileageIndex] = sol->Y[posCovered];
    sol->indexPath[bestMileageIndex] = sol->indexPath[posCovered];

    int i = posCovered;

    if (EM_USE_FAST_SOLUTION_UPDATE)
    {
        // shift elements forward of 1 position iteratively with avx until vector is too big for the amount of elements to shift (do AVX_VEC_SIZE elements per iteration)
        for (i -= AVX_VEC_SIZE; i > (int)anchorIndex; i -= AVX_VEC_SIZE)
        {
            __m256 xVec = _mm256_loadu_ps(&sol->X[i]);
            __m256 yVec = _mm256_loadu_ps(&sol->Y[i]);
            __m256i indexVec = _mm256_loadu_si256((__m256i_u*)&sol->indexPath[i]);
            _mm256_storeu_ps(&sol->X[i + 1], xVec);
            _mm256_storeu_ps(&sol->Y[i + 1], yVec);
            _mm256_storeu_si256((__m256i_u*)&sol->indexPath[i + 1], indexVec);
        }
        i += AVX_VEC_SIZE;
    }
    
    // shift elements forward one at a time
    for (i--; i > anchorIndex; i--)
    {
        sol->X[i+1] = sol->X[i];
        sol->Y[i+1] = sol->Y[i];
        sol->indexPath[i+1] = sol->indexPath[i];
    }

    i++;

    sol->X[i] = bestX;
    sol->Y[i] = bestY;
    sol->indexPath[i] = bestIndex;

    LOG(LOG_LVL_EVERYTHING, "Extra Mileage Solution Update: Node %d added to solution between nodes %d and %d", sol->indexPath[i], sol->indexPath[i-1], sol->indexPath[i+1]);
}

static void extraMileageVectorized(Solution *sol, size_t nCovered, unsigned int *rndState)
{
    // shortcuts/decluttering
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundW = inst->params.roundWeights;
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    float bestExtraMileage[AVX_VEC_SIZE];
    int bestNodes[AVX_VEC_SIZE], bestAnchors[AVX_VEC_SIZE];

    for (size_t posCovered = nCovered + 1; posCovered <= n; posCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        // Contains best mileage values
        __m256 bestExtraMileageVec = _mm256_set1_ps(INFINITY);
        // Contains the indexes of the nodes from which the best (chosen according to bestMileageVec) one will be added to the solution at the end of the iteration
        __m256i bestNodesVec = _mm256_set1_epi32(-1);
        // Contains the indexes corresponding to the edge that will be removed/ splitted to accomodate the new node
        __m256i bestAnchorsVec = _mm256_set1_epi32(-1);

        // we do this to avoid the need of checking the last elements loaded by _mm256_loadu -> exploit the "INFINITY" placed at the end of the last elements in sol.X and sol.Y
        for (size_t i = 0; i < posCovered-1; i++)
        {
            // Create vectors containig necessary data on the points attached to the edge i
            __m256 x1Vec = _mm256_broadcast_ss(&sol->X[i]), y1Vec = _mm256_broadcast_ss(&sol->Y[i]);
            __m256 x2Vec = _mm256_broadcast_ss(&sol->X[i+1]), y2Vec = _mm256_broadcast_ss(&sol->Y[i+1]);

            // Vector that contains only the cost of the current edge
            __m256 curEdgeCostVec = computeEdgeCost_VEC(x1Vec, y1Vec, x2Vec, y2Vec, ewt, roundW);

            // Vector that contains only the index of the current edge
            __m256i curEdgeID = _mm256_set1_epi32((int)i);

            // Vector that keeps track of the IDs of the best candidates for the current edge
            __m256i idsVec = _mm256_add_epi32(_mm256_set_epi32( 7, 6, 5, 4, 3, 2, 1, 0 ), _mm256_set1_epi32( posCovered ) );
            __m256i incrementVec = _mm256_set1_epi32( AVX_VEC_SIZE );

            // check for each edge which ones are the best
            for (size_t u = posCovered; u <= n; u += AVX_VEC_SIZE, idsVec = _mm256_add_epi32(idsVec, incrementVec))
            {
                __m256 curExtraMileageVec;
                {
                    __m256 xuVec = _mm256_loadu_ps(&sol->X[u]), yuVec = _mm256_loadu_ps(&sol->Y[u]);
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

        _mm256_storeu_ps(bestExtraMileage, bestExtraMileageVec);
        _mm256_storeu_si256((__m256i_u*)bestNodes, bestNodesVec);
        _mm256_storeu_si256((__m256i_u*)bestAnchors, bestAnchorsVec);

        // sort bestExtraMileage while matching the sorting moves on both bestNodes and bestAnchors
        for (size_t i = 0; i < AVX_VEC_SIZE; i++) {
            for (size_t j = i+1; j < AVX_VEC_SIZE; j++) {
                if (bestExtraMileage[i] > bestExtraMileage[j]) {
                    register float tempf;
                    swapElems(bestExtraMileage[i], bestExtraMileage[j], tempf);
                    register int tempi;
                    swapElems(bestNodes[i], bestNodes[j], tempi);
                    swapElems(bestAnchors[i], bestAnchors[j], tempi);
        } } }

        if ((inst->params.graspType != GRASP_NONE) && (rand_r(rndState) < graspThreshold) && (n - posCovered > AVX_VEC_SIZE))
        {
            switch (inst->params.graspType)
            {
            case GRASP_NONE: break;// Prevents warning in compilation

            case GRASP_ALMOSTBEST:
            {
                size_t rndIdx = 1;
                for (; rndIdx < AVX_VEC_SIZE-1; rndIdx++)
                    if (rand_r(rndState) > RAND_MAX / 2)
                        break;
                updateSolutionEM(sol, posCovered, bestNodes[rndIdx], bestAnchors[rndIdx], bestExtraMileage[rndIdx]);
                break;
            }

            case GRASP_RANDOM:
            {
                size_t node = posCovered + (size_t)rand_r(rndState) % (n+1-posCovered);
                size_t anchor = (size_t)rand_r(rndState) % (posCovered-1);
                float extraMileage = computeEdgeCost(sol->X[anchor], sol->Y[anchor], sol->X[node], sol->Y[node], ewt, roundW) +
                                     computeEdgeCost(sol->X[node], sol->Y[node], sol->X[anchor+1], sol->Y[anchor+1], ewt, roundW) -
                                     computeEdgeCost(sol->X[anchor], sol->Y[anchor], sol->X[anchor+1], sol->Y[anchor+1], ewt, roundW);
                updateSolutionEM(sol, posCovered, node, anchor, extraMileage);
                break;
            }
            }
        }
        else
            // best update
            updateSolutionEM(sol, posCovered, bestNodes[0], bestAnchors[0], bestExtraMileage[0]);
    }
}

static void extraMileageBase(Solution *sol, size_t nCovered, unsigned int *rndState, bool useCostMatrix)
{
    // shortcuts/decluttering
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundW = inst->params.roundWeights;
    int graspThreshold = (int)(inst->params.graspChance * (double)RAND_MAX);

    for (size_t posCovered = nCovered + 1; posCovered <= n; posCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        float bestMileage = INFINITY; // saves best mileage value
        size_t bestNode = 0xFFFFFFFFFFFFFFFF; // saves the index of the node that will be added at the end of the iteration

        // index of the covered nodes that represent the edge that will be splitted to integrate u at solution update
        size_t bestAnchor = 0xFFFFFFFFFFFFFFFF;

        for (size_t i = 0; i < posCovered-1; i++) // u stands for uncovered
        {
            // cost of edge already in solution [i,j]
            float currEdgeCost;
            if (useCostMatrix == 1)
                currEdgeCost = inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i + 1]];
            else
                currEdgeCost = computeEdgeCost(sol->X[i], sol->Y[i], sol->X[i+1], sol->Y[i+1], ewt, roundW);

            for (size_t u = posCovered; u < n+1; u++) // covered node i from 0 to posCovered
            {
                // sum of cost of edge [i,u] and edge [u,j]
                float altEdgeCost;
                if (useCostMatrix == 1)
                    altEdgeCost = computeEdgeCost(sol->X[u], sol->Y[u], sol->X[i], sol->Y[i], ewt, roundW) + computeEdgeCost(sol->X[u], sol->Y[u], sol->X[i+1], sol->Y[i+1], ewt, roundW);
                else
                    altEdgeCost = inst->edgeCostMat[sol->indexPath[u] * n + sol->indexPath[i]] + inst->edgeCostMat[sol->indexPath[u] * n + sol->indexPath[i+1]];

                float extraMileage = altEdgeCost - currEdgeCost;

                if (bestMileage > extraMileage)
                {
                    bestMileage = extraMileage;
                    bestNode = u;
                    bestAnchor = i;
                }
            }
            // reached this point we found the best anchors combination for the uncovered node at index u in sol.X-Y
        }
        // reached this point we found the absolute best combination of anchors and u possible for this iteration
        // we must update the solution now
        if ((inst->params.graspType != GRASP_NONE) && (rand_r(rndState) < graspThreshold) && (n - posCovered > AVX_VEC_SIZE))
        {
            switch (inst->params.graspType)
            {
            case GRASP_NONE: break; // Prevents warning in compilation
            case GRASP_ALMOSTBEST: break; // not supported in this mode
            case GRASP_RANDOM:
            {
                size_t node = posCovered + (size_t)rand_r(rndState) % (n+1-posCovered);
                size_t anchor = (size_t)rand_r(rndState) % (posCovered-1);
                float extraMileage = computeEdgeCost(sol->X[anchor], sol->Y[anchor], sol->X[node], sol->Y[node], ewt, roundW) +
                                     computeEdgeCost(sol->X[node], sol->Y[node], sol->X[anchor+1], sol->Y[anchor+1], ewt, roundW) -
                                     computeEdgeCost(sol->X[anchor], sol->Y[anchor], sol->X[anchor+1], sol->Y[anchor+1], ewt, roundW);
                updateSolutionEM(sol, posCovered, node, anchor, extraMileage);
                break;
            }
            }
        }
        else
            updateSolutionEM(sol, posCovered, bestNode, bestAnchor, bestMileage);
    }
}
