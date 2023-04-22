#include "ExtraMileage.h"
#include "EdgeCostFunctions.h"

#include "pthread.h"
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// Flag to define at compile time whether to use or not the AVX instruction to move the data
#define EM_USE_FAST_SOLUTION_UPDATE 0

// Threshold of the number of elements to move in the solution needed to use AVX to improve update speed
#define EM_FAST_SOLUTION_ELEMS_THRESHOLD 8*20


static size_t initialization(Solution *sol, enum EMInitType emType);

static void extremeCoordsPointsInit(Solution *sol);

static void farthestPointsInit(Solution *sol);

static int checkSolutionIntegrity(Solution *sol);

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2);

static inline void updateSolutionEM(Solution *sol, size_t posCovered, size_t bestMileageIndex, size_t anchorInde, float extraCost);

static void extraMileageVectorizedST(Solution *sol, size_t nCovered);

static void extraMileageST(Solution *sol, size_t nCovered);

static void extraMileageCostMatrixST(Solution *sol, size_t nCovered);


Solution ExtraMileage(Instance *inst, enum EMOptions emOpt, enum EMInitType emInitType)
{
    Solution sol = newSolution(inst);
    sol.cost = 0;

    // Set data for solution (copy coords from distance and create index path as 0,1,2,...,n-1)
    for (size_t i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
        sol.X[i] = inst->X[i];
    for (int i = 0; i < inst->nNodes; i++)
        sol.indexPath[i] = i;

    // initialize solution as specified by options
    // Set random seed for initialization
    if (inst->params.randomSeed != -1)
        srand(inst->params.randomSeed);
    else
        srand(time(NULL));
    size_t coveredNodes = initialization(&sol, emInitType);

    // apply extra mileage
    sol.execTime = applyExtraMileage(&sol, coveredNodes, emOpt);

    return sol;
}

double applyExtraMileage(Solution *sol, size_t nCovered, enum EMOptions emOpt)
{
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    // Set random seed
    if (inst->params.randomSeed != -1)
        srand(inst->params.randomSeed);
    else
        srand(time(NULL));

    // first check that the solution element at index nCovered is the same as the one at index 0
    if (sol->indexPath[0] != sol->indexPath[nCovered])
    {
        // check solution integrity when debugging
        if (LOG_LEVEL >= LOG_LVL_DEBUG || !(sol->X[n] == INFINITY && sol->Y[n] == INFINITY) || !(sol->X[n] == sol->X[0] && sol->Y[n] == sol->Y[0]))
            if (checkSolutionIntegrity(sol) != 0)
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
        extraMileageVectorizedST(sol, nCovered);
        break;
    case EM_OPTION_BASE:
        extraMileageST(sol, nCovered);
        break;
    case EM_OPTION_USE_COST_MATRIX:
        if (inst->edgeCostMat)
            extraMileageCostMatrixST(sol, nCovered);
        else
        {
            LOG(LOG_LVL_WARNING, "applyExtraMileage: edgeCostMat has not been detected/initialized. Switching to option: EM_OPTION_AVX");
            extraMileageVectorizedST(sol, nCovered);
        }
        break;
    }
    
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    double elapsed = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);
    return elapsed;
}

static size_t initialization(Solution *sol, enum EMInitType emInitType)
{
    Instance *inst = sol->instance;
    size_t coveredElems = 0;

    switch (emInitType)
    {
    case EM_INIT_RANDOM:
        // select two random nodes
        int rndIndex0 = rand() % (int)inst->nNodes, rndIndex1 = rand() % (int)inst->nNodes;
        while (abs(rndIndex1 - rndIndex0) <= 1)
            rndIndex1 = rand();

        // add the two nodes to the solution (order does not matter)
        swapElementsInSolution(sol, 0, rndIndex0);
        swapElementsInSolution(sol, 1, rndIndex1);

        // update cost
        sol->cost = computeEdgeCost(sol->X[0], sol->Y[0], sol->X[1], sol->Y[1], inst->params.edgeWeightType, inst->params.roundWeightsFlag) * 2.;

        LOG(LOG_LVL_DEBUG, "ExtraMileage-InitializationRandom: Randomly chosen edge is edge e = (%d, %d)", rndIndex0, rndIndex1);

        coveredElems = 2;
        break;
    case EM_INIT_EXTREMES:
        extremeCoordsPointsInit(sol);

        // update cost
        sol->cost = computeEdgeCost(sol->X[0], sol->Y[0], sol->X[1], sol->Y[1], inst->params.edgeWeightType, inst->params.roundWeightsFlag) * 2.;

        coveredElems = 2;
        break;
    case EM_INIT_FARTHEST_POINTS:
        farthestPointsInit(sol);

        // update cost
        sol->cost = computeEdgeCost(sol->X[0], sol->Y[0], sol->X[1], sol->Y[1], inst->params.edgeWeightType, inst->params.roundWeightsFlag) * 2.;

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

static void extremeCoordsPointsInit(Solution *sol)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
 
    __m256 maxXVec = _mm256_set1_ps(-INFINITY), maxYVec = _mm256_set1_ps(-INFINITY);
    __m256 minXVec = _mm256_set1_ps(INFINITY), minYVec = _mm256_set1_ps(INFINITY);

    __m256i maxIDVec = _mm256_set1_epi32(-1);//, maxYIDVec = _mm256_set1_epi32(-1);
    __m256i minIDVec = _mm256_set1_epi32(-1);//, minYIDVec = _mm256_set1_epi32(-1);

    __m256i idsVec = _mm256_set_epi32( 7, 6, 5, 4, 3, 2, 1, 0 ), incrementVec = _mm256_set1_epi32(AVX_VEC_SIZE);

    for (size_t i = 0; i < n; i += AVX_VEC_SIZE, idsVec = _mm256_add_epi32(idsVec, incrementVec)) // all but last iteration
    {
        if (i > n - AVX_VEC_SIZE)
        {
            idsVec = _mm256_sub_epi32(idsVec, _mm256_set1_epi32(AVX_VEC_SIZE - n % AVX_VEC_SIZE ));
            i = n - AVX_VEC_SIZE; // subtract AVX_VEC_SIZE one extra time to compensate the increment of the loop
        }

        __m256 mask, maskGT, maskLT, maskEQ;
        __m256 x = _mm256_loadu_ps(&inst->X[i]), y = _mm256_loadu_ps(&inst->Y[i]);

        // find proper mask for xMax (select only if x[j] > xMax[j] or (x[j] == xMax[j] and y[j] > yMax[j]) )
        maskEQ = _mm256_cmp_ps(maxXVec, x, _CMP_EQ_OQ);
        maskGT = _mm256_cmp_ps(y, maxYVec, _CMP_GT_OQ);
        mask = _mm256_and_ps(maskEQ, maskGT);
        maskGT = _mm256_cmp_ps(x, maxXVec, _CMP_GT_OQ);
        mask = _mm256_or_ps(mask, maskGT);

        // save maxiumum value
        maxXVec = _mm256_blendv_ps(maxXVec, x, mask);
        maxYVec = _mm256_blendv_ps(maxYVec, y, mask);
        maxIDVec = _mm256_blendv_epi8(maxIDVec, idsVec, _mm256_castps_si256(mask));

        
        // same as for xMax but for xMin
        maskEQ = _mm256_cmp_ps(minXVec, x, _CMP_EQ_OQ);
        maskLT = _mm256_cmp_ps(y, minYVec, _CMP_LT_OQ);
        mask = _mm256_and_ps(maskEQ, maskLT);
        maskLT = _mm256_cmp_ps(x, minXVec, _CMP_LT_OQ);
        mask = _mm256_or_ps(mask, maskLT);

        // save minimum values
        minXVec = _mm256_blendv_ps(minXVec, x, mask);
        minYVec = _mm256_blendv_ps(minYVec, y, mask);
        minIDVec = _mm256_blendv_epi8(minIDVec, idsVec, _mm256_castps_si256(mask));
    }

    // xMax and xMin are inside the vectors now with the corresponing index stored in the appropriate vector
    float xVecStore[AVX_VEC_SIZE];
    int IDVecStore[AVX_VEC_SIZE];

    // find and save point one
    _mm256_storeu_ps(xVecStore, maxXVec);
    _mm256_storeu_si256((__m256i*)IDVecStore, maxIDVec);
    register size_t index = 0;
    for (size_t i = 1; i < AVX_VEC_SIZE; i++)
        if (xVecStore[index] < xVecStore[i])
            index = i;
    swapElementsInSolution(sol, 0, IDVecStore[index]);

    // find and save point two
    _mm256_storeu_ps(xVecStore, minXVec);
    _mm256_storeu_si256((__m256i*)IDVecStore, minIDVec);
    index = 0;
    for (size_t i = 1; i < AVX_VEC_SIZE; i++)
        if (xVecStore[index] > xVecStore[i])
            index = i;
    swapElementsInSolution(sol, 1, IDVecStore[index]);
}

static void farthestPointsInit(Solution *sol)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    __m256 maxCostVec = _mm256_set1_ps(0), rowMaxCostVec = _mm256_set1_ps(0); // cost is always positive
    __m256i maxColIDsVec = _mm256_set1_epi32(0), maxRowIDsVec = _mm256_set1_epi32(0);
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
            __m256 costVec = computeEdgeCost_VEC(x1, y1, x2, y2, inst->params.edgeWeightType, inst->params.roundWeightsFlag);

            // check if there are costier connections in this iteration and save results
            __m256 mask = _mm256_cmp_ps(costVec, rowMaxCostVec, _CMP_GT_OQ);
            rowMaxCostVec = _mm256_blendv_ps(rowMaxCostVec, costVec, mask);
            maxColIDsVec = _mm256_blendv_epi8(maxColIDsVec, colIDsVec, _mm256_castps_si256(mask));
        }
        
        __m256 mask = _mm256_cmp_ps(rowMaxCostVec, maxCostVec, _CMP_GT_OQ);
        maxCostVec = _mm256_blendv_ps(maxCostVec, rowMaxCostVec, mask);
        maxRowIDsVec = _mm256_blendv_epi8(maxRowIDsVec, rowIDsVec, _mm256_castps_si256(mask));
        
    }

    // ####### FIND MAXIMUM COST AND CORRESPONDING INDEX USING AVX PERMUTATIONS AND COMPARISIONS

    /* permutation order
     * 4 5 6 7 3 2 1 0 
     * 2 3 1 0 6 7 5 4
     * 1 0 3 2 5 4 7 6
    */

    {   // FIRST PERMUTATION
        __m256i permuteMask = _mm256_set_epi32( 0, 1, 2, 3, 7, 6, 5, 4 );
        __m256 maxCostVecPerm = _mm256_permutevar8x32_ps(maxCostVec, permuteMask);
        __m256 mask = _mm256_cmp_ps(maxCostVecPerm, maxCostVec, _CMP_GT_OQ);
        maxCostVec = _mm256_blendv_ps(maxCostVec, maxCostVecPerm, mask);
        maxRowIDsVec = _mm256_blendv_epi8(maxRowIDsVec, _mm256_permutevar8x32_epi32(maxRowIDsVec, permuteMask), _mm256_castps_si256(mask));
        maxColIDsVec = _mm256_blendv_epi8(maxColIDsVec, _mm256_permutevar8x32_epi32(maxColIDsVec, permuteMask), _mm256_castps_si256(mask));
    }
    {   // SECOND PERMUTATION (faster instruction)
        __m256 maxCostVecPerm = _mm256_permute_ps(maxCostVec, 30); // 30 in binary 8 bits is 00011110
        __m256 mask = _mm256_cmp_ps(maxCostVecPerm, maxCostVec, _CMP_GT_OQ);
        maxCostVec = _mm256_blendv_ps(maxCostVec, maxCostVecPerm, mask);
        maxRowIDsVec = _mm256_blendv_epi8(maxRowIDsVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(maxRowIDsVec), 30)), _mm256_castps_si256(mask));
        maxColIDsVec = _mm256_blendv_epi8(maxColIDsVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(maxColIDsVec), 30)), _mm256_castps_si256(mask));
    }
    {
        // SECOND PERMUTATION (faster instruction)
        __m256 maxCostVecPerm = _mm256_permute_ps(maxCostVec, 177); // 177 in binary 8 bits is 10110001
        __m256 mask = _mm256_cmp_ps(maxCostVecPerm, maxCostVec, _CMP_GT_OQ);
        maxCostVec = _mm256_blendv_ps(maxCostVec, maxCostVecPerm, mask);
        maxRowIDsVec = _mm256_blendv_epi8(maxRowIDsVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(maxRowIDsVec), 177)), _mm256_castps_si256(mask));
        maxColIDsVec = _mm256_blendv_epi8(maxColIDsVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(maxColIDsVec), 177)), _mm256_castps_si256(mask));
    }

    // now maxCost, maxRowIDsVec, maxColIDsVec should contain only the value corresponding to the element with maximum cost so we can extract such values
    // from any position in the vector
    int maxCost = _mm_extract_ps(_mm256_castps256_ps128(maxCostVec), 0);
    LOG(LOG_LVL_DEBUG, "Extra Mileage EM_INIT_FARTHEST_POINTS: maximum cost found is %f", *(float*)&maxCost);
    int maxIndex0 = _mm_extract_epi32(_mm256_castsi256_si128(maxRowIDsVec), 0);
    int maxIndex1 = _mm_extract_epi32(_mm256_castsi256_si128(maxColIDsVec), 0);

    // initialize solution
    swapElementsInSolution(sol, 0, maxIndex0);
    swapElementsInSolution(sol, 1, maxIndex1);
}

static int checkSolutionIntegrity(Solution *sol)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    for (size_t i = 0; i < n; i++)
    {
        int index = sol->indexPath[i];
        if (index < 0 || index > n)
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: sol.indexPath[%lu] = %d which is not within the limits", i, index);
            return 1;
        }
        if (sol->X[i] != inst->X[index] || sol->Y[i] != inst->Y[index])
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: Mismatch at index %lu in solution", i);
            return 1;
        }
    }

    // everything checks out
    return 0;
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

    size_t i = posCovered;

    if (EM_USE_FAST_SOLUTION_UPDATE && (posCovered - anchorIndex > EM_FAST_SOLUTION_ELEMS_THRESHOLD))
    {
        // shift elements forward of 1 position iteratively with avx until vector is too big for the amount of elements to shift (do AVX_VEC_SIZE elements per iteration)
        for (i -= AVX_VEC_SIZE; i > anchorIndex; i -= AVX_VEC_SIZE)
        {
            __m256 xVec = _mm256_loadu_ps(&sol->X[i - AVX_VEC_SIZE]);
            __m256 yVec = _mm256_loadu_ps(&sol->Y[i - AVX_VEC_SIZE]);
            __m256i indexVec = _mm256_loadu_si256((__m256i_u*)&sol->indexPath[i - AVX_VEC_SIZE]);
            _mm256_storeu_ps(&sol->X[i - AVX_VEC_SIZE + 1], xVec);
            _mm256_storeu_ps(&sol->Y[i - AVX_VEC_SIZE + 1], yVec);
            _mm256_storeu_si256((__m256i_u*)&sol->indexPath[i - AVX_VEC_SIZE + 1], indexVec);
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

static void extraMileageVectorizedST(Solution *sol, size_t nCovered)
{
    // shortcuts/decluttering
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    enum edgeWeightType ewt = inst->params.edgeWeightType;
    int roundW = inst->params.roundWeightsFlag;

    float bestCostVecStore[8] = { 0 };
    size_t bestCostSortedIndexes[8] = { 0 };

    for (size_t posCovered = nCovered + 1; posCovered <= n; posCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        // Contains best mileage values
        __m256 bestExtraMileageVec = _mm256_set1_ps(INFINITY);
        // Contains the indexes of the nodes from which the best (chosen according to bestMileageVec) one will be added to the solution at the end of the iteration
        __m256i bestExtraMileageNodeIDVec = _mm256_set1_epi32(-1);
        // Contains the indexes corresponding to the edge that will be removed/ splitted to accomodate the new node
        __m256i bestExtraMileageAnchorIDVec = _mm256_set1_epi32(-1);

        // approach is opposite compared to the non-vectorized approach: check for each edge (i) in the tour all uncovered nodes
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
                //{
                    __m256 xuVec = _mm256_loadu_ps(&sol->X[u]), yuVec = _mm256_loadu_ps(&sol->Y[u]);
                    __m256 altEdge1CostVec = computeEdgeCost_VEC(xuVec, yuVec, x1Vec, y1Vec, ewt, roundW);
                    __m256 altEdge2CostVec = computeEdgeCost_VEC(xuVec, yuVec, x2Vec, y2Vec, ewt, roundW);
                    __m256 altEdgeCostVec = _mm256_add_ps(altEdge1CostVec, altEdge2CostVec);
                    curExtraMileageVec = _mm256_sub_ps(altEdgeCostVec, curEdgeCostVec);
                //}

                // Compare curExtraMileageCostVec with bestExtraMileageVec
                __m256 cmpMask = _mm256_cmp_ps(curExtraMileageVec, bestExtraMileageVec, _CMP_LT_OQ);

                // Set new best according to comparison result
                bestExtraMileageVec = _mm256_blendv_ps(bestExtraMileageVec, curExtraMileageVec, cmpMask);
                bestExtraMileageAnchorIDVec = _mm256_blendv_epi8(bestExtraMileageAnchorIDVec, curEdgeID, _mm256_castps_si256(cmpMask));
                bestExtraMileageNodeIDVec = _mm256_blendv_epi8(bestExtraMileageNodeIDVec, idsVec, _mm256_castps_si256(cmpMask));
            }
        }
        // at this point we must select the best canditate(the one with minimum cost) in bestExtraMileageVec

        _mm256_storeu_ps(bestCostVecStore, bestExtraMileageVec);
        for (size_t sortedElemsCount = 0; sortedElemsCount < AVX_VEC_SIZE; sortedElemsCount++)
        {
            float min = INFINITY;
            for (size_t i = 0; i < AVX_VEC_SIZE; i++)
            {
                int flag = 1;
                for (size_t j = 0; j < sortedElemsCount; j++)
                    if (bestCostSortedIndexes[j] == i)
                        flag = 0;

                if (flag && min > bestCostVecStore[i])
                {
                    min = bestCostVecStore[i];
                    bestCostSortedIndexes[sortedElemsCount] = i;
                }
            }
        }
        


        // we do this by comparing the vector with its permutation multiple times, until all the elements in the vec are equal to the best one

        { // FIRST PERMUTATION
            __m256i permuteMask = _mm256_set_epi32( 0, 1, 2, 3, 7, 6, 5, 4 );
            __m256 permutedBestExtraMilageVec = _mm256_permutevar8x32_ps(bestExtraMileageVec, permuteMask);
            __m256 mask = _mm256_cmp_ps(permutedBestExtraMilageVec, bestExtraMileageVec, _CMP_LT_OQ);
            bestExtraMileageVec = _mm256_blendv_ps(bestExtraMileageVec, permutedBestExtraMilageVec, mask);
            bestExtraMileageNodeIDVec = _mm256_blendv_epi8(bestExtraMileageNodeIDVec, _mm256_permutevar8x32_epi32(bestExtraMileageNodeIDVec, permuteMask), _mm256_castps_si256(mask));
            bestExtraMileageAnchorIDVec = _mm256_blendv_epi8(bestExtraMileageAnchorIDVec, _mm256_permutevar8x32_epi32(bestExtraMileageAnchorIDVec, permuteMask), _mm256_castps_si256(mask));
        }
        { // SECOND PERMUTATION
            __m256 permutedBestExtraMilageVec = _mm256_permute_ps(bestExtraMileageVec, 30); // 30 in binary 8 bits is 00011110
            __m256 mask = _mm256_cmp_ps(permutedBestExtraMilageVec, bestExtraMileageVec, _CMP_LT_OQ);
            bestExtraMileageVec = _mm256_blendv_ps(bestExtraMileageVec, permutedBestExtraMilageVec, mask);
            bestExtraMileageNodeIDVec = _mm256_blendv_epi8(bestExtraMileageNodeIDVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(bestExtraMileageNodeIDVec), 30)), _mm256_castps_si256(mask));
            bestExtraMileageAnchorIDVec = _mm256_blendv_epi8(bestExtraMileageAnchorIDVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(bestExtraMileageAnchorIDVec), 30)), _mm256_castps_si256(mask));
        }
        { // THIRD PERMUTATION
            __m256 permutedBestExtraMilageVec = _mm256_permute_ps(bestExtraMileageVec, 177); // 177 in binary 8 bits is 10110001
            __m256 mask = _mm256_cmp_ps(permutedBestExtraMilageVec, bestExtraMileageVec, _CMP_LT_OQ);
            bestExtraMileageVec = _mm256_blendv_ps(bestExtraMileageVec, permutedBestExtraMilageVec, mask);
            bestExtraMileageNodeIDVec = _mm256_blendv_epi8(bestExtraMileageNodeIDVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(bestExtraMileageNodeIDVec), 177)), _mm256_castps_si256(mask));
            bestExtraMileageAnchorIDVec = _mm256_blendv_epi8(bestExtraMileageAnchorIDVec, _mm256_castps_si256(_mm256_permute_ps(_mm256_castsi256_ps(bestExtraMileageAnchorIDVec), 177)), _mm256_castps_si256(mask));
        }

        int bestCostInt = _mm_extract_epi32(_mm256_castsi256_si128(_mm256_castps_si256(bestExtraMileageVec)), 0);
        float bestCost = *(float*)&bestCostInt;
        size_t bestExtraMileageNodeID = (size_t)_mm_extract_epi32(_mm256_castsi256_si128(bestExtraMileageNodeIDVec), 0);
        size_t bestExtraMileageAnchorID = (size_t)_mm_extract_epi32(_mm256_castsi256_si128(bestExtraMileageAnchorIDVec), 0);

        // reached this point we found the absolute best combination of anchors and u possible for this iteration
        // we must update the solution now
        updateSolutionEM(sol, posCovered, bestExtraMileageNodeID, bestExtraMileageAnchorID, bestCost);
    }
}

static void extraMileageST(Solution *sol, size_t nCovered)
{
    // shortcuts/decluttering
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    enum edgeWeightType ewt = inst->params.edgeWeightType;
    int roundW = inst->params.roundWeightsFlag;

    for (size_t posCovered = nCovered + 1; posCovered <= n; posCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        float bestMileage = INFINITY; // saves best mileage value
        size_t bestMileageNodeID = 0xFFFFFFFFFFFFFFFF; // saves the index of the node that will be added at the end of the iteration

        // index of the covered nodes that represent the edge that will be splitted to integrate u at solution update
        size_t bestMileageAnchor = 0xFFFFFFFFFFFFFFFF;

        for (size_t i = 0; i < posCovered-1; i++) // u stands for uncovered
        {
            // cost of edge already in solution [i,j]
            float currEdgeCost = computeEdgeCost(sol->X[i], sol->Y[i], sol->X[i+1], sol->Y[i+1], ewt, roundW);

            for (size_t u = posCovered; u < n+1; u++) // covered node i from 0 to posCovered
            {
                // cost of edge [i,u]
                float cost1 = computeEdgeCost(sol->X[u], sol->Y[u], sol->X[i], sol->Y[i], ewt, roundW);
                // sum of cost of edge [i,u] and edge [u,j]
                float altEdgeCost = cost1 + computeEdgeCost(sol->X[u], sol->Y[u], sol->X[i+1], sol->Y[i+1], ewt, roundW);

                float extraCost = altEdgeCost - currEdgeCost;

                if (bestMileage > extraCost)
                {
                    bestMileage = extraCost;
                    bestMileageNodeID = u;
                    bestMileageAnchor = i;
                }
            }
            // reached this point we found the best anchors combination for the uncovered node at index u in sol.X-Y
        }
        // reached this point we found the absolute best combination of anchors and u possible for this iteration
        // we must update the solution now
        updateSolutionEM(sol, posCovered, bestMileageNodeID, bestMileageAnchor, bestMileage);
    }
}


static void extraMileageCostMatrixST(Solution *sol, size_t nCovered)
{
    // shortcuts/decluttering
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    float *edgeCostMat = inst->edgeCostMat;

    for (size_t posCovered = nCovered + 1; posCovered <= n; posCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        float bestMileage = INFINITY; // saves best mileage value
        size_t bestMileageNodeID = 0xFFFFFFFFFFFFFFFF; // saves the index of the node that will be added at the end of the iteration

        // index of the covered nodes that represent the edge that will be splitted to integrate u at solution update
        size_t bestMileageAnchor = 0xFFFFFFFFFFFFFFFF;

        for (size_t i = 0; i < posCovered-1; i++) // u stands for uncovered
        {
            // cost of edge already in solution [i,j]
            float currEdgeCost = edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i + 1]];

            for (size_t u = posCovered; u <= n; u++) // covered node i from 0 to posCovered
            {
                // cost of edge [i,u]
                float cost1 = edgeCostMat[sol->indexPath[u] * n + sol->indexPath[i]];
                // sum of cost of edge [i,u] and edge [u,j]
                float altEdgeCost = cost1 + edgeCostMat[sol->indexPath[u] * n + sol->indexPath[i+1]];

                float extraCost = altEdgeCost - currEdgeCost;

                if (bestMileage > extraCost)
                {
                    bestMileage = extraCost;
                    bestMileageNodeID = u;
                    bestMileageAnchor = i;
                }
            }
            // reached this point we found the best anchors combination for the uncovered node at index u in sol.X-Y
        }
        // reached this point we found the absolute best combination of anchors and u possible for this iteration
        // we must update the solution now
        updateSolutionEM(sol, posCovered, bestMileageNodeID, bestMileageAnchor, bestMileage);
    }
}

