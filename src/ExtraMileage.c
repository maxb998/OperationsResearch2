#include "ExtraMileage.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"

#include "pthread.h"
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

#define EM_USE_FAST_SOLUTION_UPDATE 0
#define EM_FAST_SOLUTION_ELEMS_THRESHOLD 500


static size_t initialization(Solution *sol, enum EMInitType emType);

static void extremeCoordsPointsInit(Solution *sol);

static void farthestPointsInit(Solution *sol);

static inline void swapSolutions(Solution *s1, Solution *s2);

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2);

static inline void updateSolutionExtraMileage(Solution *sol, size_t posCovered, size_t bestMileageIndex, size_t anchorIndex);

static void extraMileageST(Solution *sol, size_t nCovered);

static void extraMileageVectorizedST(Solution *sol, enum EMInitType emType);



Solution ExtraMileage(Instance *inst, enum EMInitType emType)
{
    Solution sol = newSolution(inst);

    // Set data for solution (copy coords from distance and create index path as 0,1,2,...,n-1)
    for (size_t i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
        sol.X[i] = inst->X[i];
    for (int i = 0; i < inst->nNodes; i++)
        sol.indexPath[i] = i;

    size_t coveredNodes = initialization(&sol, emType);

    extraMileageST(&sol, coveredNodes);

    return sol;
}

void runExtraMileageOnce(Solution *sol, size_t nCovered)
{
    size_t n = sol->instance->nNodes;

    for (size_t u = nCovered+1; u <= n; u++)
    {
        for (size_t i = 0; i <= nCovered+1; i++)
        {
            /* code */
        }
        
    }
    
}

static size_t initialization(Solution *sol, enum EMInitType emType)
{
    Instance *inst = sol->instance;
    size_t coveredElems = 0;

    switch (emType)
    {
    case EM_INIT_RANDOM:
        // select two random nodes
        int rndIndex0 = rand() % (int)inst->nNodes, rndIndex1 = rand() % (int)inst->nNodes;
        while (abs(rndIndex1 - rndIndex0) <= 1)
            rndIndex1 = rand();

        // add the two nodes to the solution (order does not matter)
        swapElementsInSolution(sol, 0, rndIndex0);
        swapElementsInSolution(sol, 1, rndIndex1);

        coveredElems = 2;
        break;
    case EM_INIT_EXTREMES:
        extremeCoordsPointsInit(sol);
        coveredElems = 2;
        break;
    case EM_INIT_FARTHEST_POINTS:
        farthestPointsInit(sol);
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
    //sol->X[0] = inst->X[IDVecStore[index]];
    //sol->Y[0] = inst->Y[IDVecStore[index]];
    //sol->indexPath[0] = IDVecStore[index];

    // find and save point two
    _mm256_storeu_ps(xVecStore, minXVec);
    _mm256_storeu_si256((__m256i*)IDVecStore, minIDVec);
    index = 0;
    for (size_t i = 1; i < AVX_VEC_SIZE; i++)
        if (xVecStore[index] > xVecStore[i])
            index = i;
    swapElementsInSolution(sol, 1, IDVecStore[index]);
    //sol->X[1] = inst->X[IDVecStore[index]];
    //sol->Y[1] = inst->Y[IDVecStore[index]];
    //sol->indexPath[1] = IDVecStore[index];
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
            __m256 costVec = computeEdgeCost_VEC(x1, y1, x2, y2, inst->params.edgeWeightType, inst->params.roundWeights);

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
    //_MM_EXTRACT_FLOAT();
    int maxCost = _mm_extract_ps(_mm256_castps256_ps128(maxCostVec), 0);
    LOG(LOG_LVL_DEBUG, "Extra Mileage EM_INIT_FARTHEST_POINTS: maximum cost found is %f", *(float*)&maxCost);
    int maxIndex0 = _mm_extract_epi32(_mm256_castsi256_si128(maxRowIDsVec), 0);
    int maxIndex1 = _mm_extract_epi32(_mm256_castsi256_si128(maxColIDsVec), 0);

    // initialize solution
    swapElementsInSolution(sol, 0, maxIndex0);
    swapElementsInSolution(sol, 1, maxIndex1);
    /*sol->indexPath[0] = maxIndex0;
    sol->indexPath[1] = maxIndex1;
    sol->X[0] = inst->X[maxIndex0];
    sol->X[1] = inst->X[maxIndex1];
    sol->Y[0] = inst->Y[maxIndex0];
    sol->Y[1] = inst->Y[maxIndex1];// */
}

static inline void swapSolutions(Solution *s1, Solution *s2)
{
    register float swapCost;
    swapElems(s1->bestCost, s2->bestCost, swapCost);

    register float *swapf;
    swapElems(s1->X, s2->X, swapf);
    swapElems(s1->Y, s2->Y, swapf);

    register int *swapi;
    swapElems(s1->indexPath, s2->indexPath, swapi);
}

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2)
{
    register float tempf;
    swapElems(sol->X[pos1], sol->X[pos2], tempf);
    swapElems(sol->Y[pos1], sol->Y[pos2], tempf);
    register int tempi;
    swapElems(sol->indexPath[pos1], sol->indexPath[pos2], tempi);
}

static inline void updateSolutionExtraMileage(Solution *sol, size_t posCovered, size_t bestMileageIndex, size_t anchorIndex)
{
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

    LOG(LOG_LVL_EVERYTHING, "Extra Mileage Solution Update: Node %d added to solution between nodes %d and %d",
                                                    sol->indexPath[i], sol->indexPath[i-1], sol->indexPath[i+1]);
}

static void extraMileageST(Solution *sol, size_t nCovered)
{
    // shortcuts/decluttering
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;
    enum edgeWeightType ewt = inst->params.edgeWeightType;
    int roundW = inst->params.roundWeights;

    for (size_t posCovered = nCovered + 1; posCovered <= n; posCovered++) // until there are uncored nodes (each iteration adds one to posCovered)
    {
        float bestMileage = INFINITY; // saves best mileage value
        size_t bestMileageNodeID = 0xFFFFFFFFFFFFFFFF; // saves the index of the node that will be added at the end of the iteration

        // index of the covered nodes that represent the edge that will be splitted to integrate u at solution update
        size_t bestMileageAnchor = 0xFFFFFFFFFFFFFFFF;

        for (size_t u = posCovered; u < n+1; u++) // u stands for uncovered
        {
            for (size_t i = 0; i < posCovered-1; i++) // covered node i from 0 to posCovered
            {
                // cost of edge already in solution [i,j]
                float currEdgeCost = computeEdgeCost(sol->X[i], sol->Y[i], sol->X[i+1], sol->Y[i+1], ewt, roundW);
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
        updateSolutionExtraMileage(sol, posCovered, bestMileageNodeID, bestMileageAnchor);
    }
    //sol->X[n] = sol->X[0];
    //sol->Y[n] = sol->Y[0];
    //sol->indexPath[n] = sol->indexPath[0];
}


