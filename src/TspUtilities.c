#include "TspUtilities.h"
#include "EdgeCostFunctions.h"


// Check the correctness of the cost of the solution stored in Solution sol.
static void checkCost(Solution *sol);

void checkSolution(Solution *sol)
{
    Instance *inst = sol->instance;

    char * coveredNodes = calloc(sol->instance->nNodes, sizeof(char));

    // First and last node must be equal (the circuit is closed)
    if (sol->indexPath[0] != sol->indexPath[sol->instance->nNodes]) 
        throwError(sol->instance, sol, "SolutionCheck: first and last node in sol.indexPath should coincide, but they do not");
    // also check for sol.X and sol.Y
    if (sol->X[0] != sol->X[inst->nNodes])
        throwError(inst, sol, "SolutionCheck: first and last node in sol.X should coincide, but they do not");
    if (sol->Y[0] != sol->Y[inst->nNodes])
        throwError(inst, sol, "SolutionCheck: first and last node in sol.Y should coincide, but they do not");

    LOG(LOG_LVL_EVERYTHING, "SolutionCheck: first and last node in sol.indexPath coincide.");

    // Populate uncoveredNodes array, here we check if a node is repeated along the path
    for (int i = 0; i < inst->nNodes; i++)
    {
        int currentNode = sol->indexPath[i];

        if(coveredNodes[currentNode] == 1)
            throwError(inst, sol, "SolutionCheck: node %d repeated in the solution. Loop iteration %d", currentNode, i);
        else
            coveredNodes[currentNode] = 1;
    }
    LOG(LOG_LVL_EVERYTHING, "SolutionCheck: all nodes in the path are unique.");

    // Check that all the nodes are covered in the path
    for (int i = 0; i < inst->nNodes; i++)
        if(coveredNodes[i] == 0)
            throwError(inst, sol, "SolutionCheck: node %d is not in the path", i);
    free(coveredNodes);
    LOG(LOG_LVL_EVERYTHING, "SolutionCheck: all the nodes are present in the path");

    // Check sol.X and sol.Y if they correspond correctly to inst.X and inst.Y given the indexes in sol.indexPath
    for (size_t i = 0; i < inst->nNodes; i++)
    {
        if (sol->X[i] != inst->X[sol->indexPath[i]])
            throwError(inst, sol, "SolutionCheck: sol.X[sol.indexPath[%ld]] = %.3e and does not correspond with inst.X[%ld] = %.3e", i, inst->X[sol->indexPath[i]], i, sol->X[i]);
        if (sol->Y[i] != inst->Y[sol->indexPath[i]])
            throwError(inst, sol, "SolutionCheck: sol.Y[sol.indexPath[%ld]] = %.3e and does not correspond with inst.Y[%ld] = %.3e", i, inst->Y[sol->indexPath[i]], i, sol->Y[i]);
    }

    checkCost(sol);

    LOG(LOG_LVL_DEBUG, "SolutionCheck: solution is coherent and feasible");
}


static void checkCost(Solution *sol)
{
    double recomputedCost = computeSolutionCostVectorizedDouble(sol);

    if (recomputedCost != sol->cost)
        throwError(sol->instance, sol, "CheckCost: Error in the computation of the pathCost. Recomputed Cost: %lf Cost in Solution: %lf", recomputedCost, sol->cost);
    
    LOG(LOG_LVL_EVERYTHING, "CheckCost: Cost of solution is correct");
}




float computeSolutionCostVectorizedFloat(Solution *sol)
{
    enum edgeWeightType ewt = sol->instance->params.edgeWeightType;
    int roundFlag = sol->instance->params.roundWeightsFlag;

    register __m256 costVec = _mm256_setzero_ps();
    size_t i = 0;
    while (i < sol->instance->nNodes - AVX_VEC_SIZE)
    {
        register __m256 x1, x2, y1, y2;
        x1 = _mm256_loadu_ps(&sol->X[i]), y1 = _mm256_loadu_ps(&sol->Y[i]);
        x2 = _mm256_loadu_ps(&sol->X[i+1]), y2 = _mm256_loadu_ps(&sol->Y[i+1]);

        costVec = _mm256_add_ps(costVec, computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag));

        i += AVX_VEC_SIZE;
    }

    // here we do the last iteration. since all "extra" elements at the end of X and Y are NaNs they can cause issues on sum. so we use a maskload for the last vector
    register __m256 x1, x2, y1, y2;
    x1 = _mm256_loadu_ps(&sol->X[i]); y1 = _mm256_loadu_ps(&sol->Y[i]);
    x2 = _mm256_loadu_ps(&sol->X[i+1]); y2 = _mm256_loadu_ps(&sol->Y[i+1]);

    register __m256 lastDist = computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag);

    // now we sustitute the infinity and NaN in the vector with zeroes
    register __m256 infinityVec = _mm256_set1_ps(INFINITY);
    register __m256 mask = _mm256_cmp_ps(lastDist, infinityVec, _CMP_LT_OQ);
    costVec = _mm256_add_ps(costVec, _mm256_blendv_ps(mask, computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag), mask));


    float costVecStore[8];
    _mm256_storeu_ps(costVecStore, costVec);
    
    double totalCost = 0.0;
    for (size_t i = 0; i < AVX_VEC_SIZE; i++)
        totalCost += costVecStore[i];
    
    return totalCost;
}

double computeSolutionCostVectorizedDouble(Solution *sol)
{
    size_t n = sol->instance->nNodes;
    enum edgeWeightType ewt = sol->instance->params.edgeWeightType;
    int roundFlag = sol->instance->params.roundWeightsFlag;

    __m256d costVec = _mm256_setzero_pd();
    size_t i = 0;
    while (i < n) //- AVX_VEC_SIZE)
    {
        __m256 x1, x2, y1, y2;

        x1 = _mm256_loadu_ps(&sol->X[i]);
        y1 = _mm256_loadu_ps(&sol->Y[i]);
        x2 = _mm256_loadu_ps(&sol->X[i+1]);
        y2 = _mm256_loadu_ps(&sol->Y[i+1]);

        if ((i > n - AVX_VEC_SIZE) && (n % AVX_VEC_SIZE != 0))
        {
            int loadMask[AVX_VEC_SIZE] = { 0 };
            for (size_t j = 0; j < n % AVX_VEC_SIZE; j++)
                loadMask[j] = -1;
            __m256i mask = _mm256_loadu_si256((__m256i_u*)loadMask);
            
            x1 = _mm256_maskload_ps(&sol->X[i], mask);
            y1 = _mm256_maskload_ps(&sol->Y[i], mask);
            x2 = _mm256_maskload_ps(&sol->X[i+1], mask);
            y2 = _mm256_maskload_ps(&sol->Y[i+1], mask);
        }

        __m256 costVecFloat = computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag);

        // convert vector of floats into 2 vectors of doubles and add them to the total cost
        // first half of the vector
        __m128 partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 0);
        __m256d partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
        costVec = _mm256_add_pd(costVec, partOfCostVecDouble);
        // second half of the vector
        partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 1);
        partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
        costVec = _mm256_add_pd(costVec, partOfCostVecDouble);

        i += AVX_VEC_SIZE;
    }
    /*
    // here we do the last iteration. since all "extra" elements at the end of X and Y are NaNs they can cause issues on sum. so we use a maskload for the last vector
    register __m256 x1, x2, y1, y2;
    x1 = _mm256_loadu_ps(&sol->X[i]); y1 = _mm256_loadu_ps(&sol->Y[i]);
    x2 = _mm256_loadu_ps(&sol->X[i+1]); y2 = _mm256_loadu_ps(&sol->Y[i+1]);

    register __m256 lastDist = computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag);

    // now we sustitute the infinities in the vector with zeroes
    register __m256 infinityVec = _mm256_set1_ps(INFINITY);
    register __m256 mask = _mm256_cmp_ps(lastDist, infinityVec, _CMP_LT_OQ);
    register __m256 costVecFloat = _mm256_blendv_ps(mask, computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag), mask);

    // convert vector of floats into 2 vectors of doubles and add them to the total cost
    // first half of the vector
    register __m128 partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 0);
    register __m256d partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
    costVec = _mm256_add_pd(costVec, partOfCostVecDouble);
    // second half of the vector
    partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 1);
    partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
    costVec = _mm256_add_pd(costVec, partOfCostVecDouble);*/

    double vecStore[4];
    _mm256_storeu_pd(vecStore, costVec);
    
    double totalCost = 0.0;
    for (size_t i = 0; i < 4; i++)
        totalCost += vecStore[i];
    
    return totalCost;
}

double computeSolutionCost(Solution *sol)
{
    enum edgeWeightType ewt = sol->instance->params.edgeWeightType;
    int roundFlag = sol->instance->params.roundWeightsFlag;

    double cost = 0.0;
    for (size_t i = 0; i < sol->instance->nNodes; i++)
        cost += (double)computeEdgeCost(sol->X[i], sol->Y[i], sol->X[i+1], sol->Y[i+1], ewt, roundFlag);
    
    return cost;
}