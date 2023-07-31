#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

#define LOG_INTERVAL 5


typedef struct
{
    float costOffset;
    int edge0;
    int edge1;
} _2optMoveData;

static bool printPerformanceLog = false;

void set2OptPerformanceBenchmarkLog(bool val)
{
    printPerformanceLog = val;
}


static inline void updateSolution(Solution *sol, _2optMoveData bestFix, float *X, float *Y);

static inline _2optMoveData _2optBestFixBase(Solution *sol, bool useCostMatrix);

static inline _2optMoveData _2OptBestFixAVX(Solution *sol, float *X, float *Y);

void apply2OptBestFix(Solution *sol)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    int iterNum = 0;

    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    /*if (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
    {
        if (inst->edgeCostMat == NULL)
        {
            LOG(LOG_LVL_WARNING, "2Opt: edgeCostMat has not been detected. Computing it now");
            double costMatTime = computeCostMatrix(inst);
            LOG(LOG_LVL_WARNING, "2Opt: edgeCostMat computation took: %lf seconds", costMatTime);
        }
    }*/

    // check integrity if debugging
    if (inst->params.logLevel >= LOG_LVL_DEBUG)
    {
        // always check solution
        if (!checkSolution(sol))
            throwError(inst, sol, "apply2OptBestFix: Input solution is not valid");
    }

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double printTimeSec = cvtTimespec2Double(timeStruct);

    float *X = NULL, *Y = NULL;

    sol->indexPath[n] = sol->indexPath[0];

    if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
    {
        X = malloc((n + AVX_VEC_SIZE) * 2 * sizeof(float));
        if (X == NULL)
            throwError(inst, sol, "_2optBestFixBase : Failed to allocate memory");
        Y = &X[n + AVX_VEC_SIZE];
        for (int i = 0; i < n; i++)
        {
            X[i] = inst->X[sol->indexPath[i]];
            Y[i] = inst->Y[sol->indexPath[i]];
        }
        X[n] = X[0];
        Y[n] = Y[0];
        for (int i = n + 1; i < n + AVX_VEC_SIZE; i++)
        {
            X[i] = INFINITY;
            Y[i] = INFINITY;
        }
    }

    bool notFinishedFlag = true;
    while (notFinishedFlag) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        // setup local values to avoid calling mutex too often
        _2optMoveData bestFix;

        switch (COMPUTATION_TYPE)
        {
        case COMPUTE_OPTION_AVX:
            bestFix = _2OptBestFixAVX(sol, X, Y);
            break;

        case COMPUTE_OPTION_BASE:
            bestFix = _2optBestFixBase(sol, false);
            break;

        case COMPUTE_OPTION_USE_COST_MATRIX:
            bestFix = _2optBestFixBase(sol, true);
            break;
        }

        if (bestFix.costOffset < -EPSILON)
            updateSolution(sol, bestFix, X, Y);
        else
            notFinishedFlag = false;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        double currentTime = cvtTimespec2Double(timeStruct);
        if (currentTime - printTimeSec > LOG_INTERVAL)
        {
            LOG(LOG_LVL_LOG, "Solution optimization in progress: cost is %lf at iteration %4lu with last optimization of %f", sol->cost, iterNum, -bestFix.costOffset);
            printTimeSec = currentTime;
        }
        iterNum++;
    }

    free(X);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double elapsed = cvtTimespec2Double(timeStruct) - startTime;
    if (printPerformanceLog)
    {
        LOG(LOG_LVL_NOTICE, "Total number of iterations: %lu", iterNum);
        LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterNum/elapsed);
    }

    sol->execTime += elapsed;
}

static inline void updateSolution(Solution *sol, _2optMoveData bestFix, float *X, float *Y)
{
    // update cost
    sol->cost += (double)bestFix.costOffset;

    /*
     *      bestSolIDs = { 0 1 2 3 4 5 6 7 8 9 }
     * if bestSolution = { 2 5 8 7 6 9 4 1 0 3 } and bestOffsetEdges = { 3 8 }
     *          -> means edges to swapped are (7,6) and (0,3) with (7,0) and (0,6)
     * at the end of the swap, the new bestSolution will be { 2 5 8 7 0 1 4 9 6 3 }     ^         ^
     *     old bestSolution = { 2 5 8 7 6 9 4 1 0 3 }   with original indexes = { 0 1 2 3 4 5 6 7 8 9 }
     *     new bestSolution = { 2 5 8 7 0 1 4 9 6 3 }   with original indexes = { 0 1 2 3 8 7 6 5 4 9 }
     *
     * Which means that we must invert the elements of bestSolution from index 3(not inlcuded) to index 8(included)
     */

    int smallID = bestFix.edge0 + 1, bigID = bestFix.edge1;

    static int updateCount = 0; // can generate issues with vns
    updateCount++;

    LOG(LOG_LVL_EVERYTHING, "2Opt: [%d] Updating solution by switching edge (%d,%d) with edge (%d,%d) improving cost by %f", updateCount,
        sol->indexPath[bestFix.edge0], sol->indexPath[bestFix.edge0 + 1],
        sol->indexPath[bestFix.edge1], sol->indexPath[bestFix.edge1 + 1], bestFix.costOffset);


    while (smallID < bigID)
    {
        register int tempInt;
        swapElems(sol->indexPath[smallID], sol->indexPath[bigID], tempInt);
        if (X)
        {
            register float tempFloat;
            swapElems(X[smallID], X[bigID], tempFloat);
            swapElems(Y[smallID], Y[bigID], tempFloat);
        }

        smallID++;
        bigID--;
    }
}

static inline _2optMoveData _2optBestFixBase(Solution *sol, bool useCostMatrix)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType;
    bool roundW = inst->params.roundWeights;

    _2optMoveData bestFix = { .costOffset=0 };

    for (_2optMoveData currFix = { .edge0=0 }; currFix.edge0 < n - 1; currFix.edge0++) // check for one edge at a time every other edge(except already checked)
    {
        float partSolEdgeWgt;
        if (useCostMatrix)
            partSolEdgeWgt = inst->edgeCostMat[sol->indexPath[currFix.edge0] * n + sol->indexPath[currFix.edge0 + 1]];
        else
            partSolEdgeWgt = computeEdgeCost(inst->X[sol->indexPath[currFix.edge0]], inst->Y[sol->indexPath[currFix.edge0]], inst->X[sol->indexPath[currFix.edge0 + 1]], inst->Y[sol->indexPath[currFix.edge0 + 1]], ewt, roundW);

        for (currFix.edge1 = 2 + currFix.edge0; (currFix.edge1 < n - 1) || ((currFix.edge1 < n) && (currFix.edge0 > 0)); currFix.edge1++)
        {
            float solEdgeWgt;
            if (useCostMatrix)
                solEdgeWgt = partSolEdgeWgt + inst->edgeCostMat[sol->indexPath[currFix.edge1] * n + sol->indexPath[currFix.edge1 + 1]];
            else
                solEdgeWgt = partSolEdgeWgt + computeEdgeCost(inst->X[sol->indexPath[currFix.edge1]], inst->Y[sol->indexPath[currFix.edge1]], inst->X[sol->indexPath[currFix.edge1 + 1]], inst->Y[sol->indexPath[currFix.edge1 + 1]], ewt, roundW);

            // check the combined weight other combination of edges
            float altEdgeWgt;
            if (useCostMatrix)
                altEdgeWgt = inst->edgeCostMat[sol->indexPath[currFix.edge0] * n + sol->indexPath[currFix.edge1]] + inst->edgeCostMat[sol->indexPath[currFix.edge1 + 1] * n + sol->indexPath[currFix.edge0 + 1]];
            else
                altEdgeWgt = computeEdgeCost(inst->X[sol->indexPath[currFix.edge0]], inst->Y[sol->indexPath[currFix.edge0]], inst->X[sol->indexPath[currFix.edge1]], inst->Y[sol->indexPath[currFix.edge1]], ewt, roundW) + 
                             computeEdgeCost(inst->X[sol->indexPath[currFix.edge0 + 1]], inst->Y[sol->indexPath[currFix.edge0 + 1]], inst->X[sol->indexPath[currFix.edge1 + 1]], inst->Y[sol->indexPath[currFix.edge1 + 1]], ewt, roundW);

            currFix.costOffset = altEdgeWgt - solEdgeWgt;
            // update local best if current one is better
            if (bestFix.costOffset > currFix.costOffset)
                bestFix = currFix;
        }
    }

    return bestFix;
}

static inline _2optMoveData _2OptBestFixAVX(Solution *sol, float *X, float *Y)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType;
    bool roundW = inst->params.roundWeights;

    _2optMoveData bestFix = { .costOffset=0 };

    for (int edge0 = 0; edge0 < n - 1; edge0++) // check for one edge at a time every other edge(except already checked)
    {

        __m256 x1 = _mm256_broadcast_ss(&X[edge0]), y1 = _mm256_broadcast_ss(&Y[edge0]);
        __m256 x2 = _mm256_broadcast_ss(&X[edge0 + 1]), y2 = _mm256_broadcast_ss(&Y[edge0 + 1]);
        __m256 partialSolEdgeWgt = computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundW);
        __m256 bestOffsetVec = _mm256_set1_ps(INFINITY);

        __m256i idsVec = _mm256_add_epi32((_mm256_set_epi32(9, 8, 7, 6, 5, 4, 3, 2)), _mm256_set1_epi32(edge0));
        __m256i increment = _mm256_set1_epi32(AVX_VEC_SIZE);
        __m256i bestIDsVec = _mm256_set1_epi32(-1);

        for (int i = 2 + edge0; (i < n - 1) || ((i < n) && (edge0 > 0)); i += AVX_VEC_SIZE)
        {
            __m256 altEdgeWgt, solEdgeWgt;
            { // scope "force" a thing compiler should do automatically -> x3,y3,x4,y4 destroyed as soon as we don't need them anymore
                __m256 x3 = _mm256_loadu_ps(&X[i]), y3 = _mm256_loadu_ps(&Y[i]);
                __m256 x4 = _mm256_loadu_ps(&X[i + 1]), y4 = _mm256_loadu_ps(&Y[i + 1]);
                solEdgeWgt = _mm256_add_ps(partialSolEdgeWgt, computeEdgeCost_VEC(x3, y3, x4, y4, ewt, roundW));

                altEdgeWgt = _mm256_add_ps(computeEdgeCost_VEC(x1, y1, x3, y3, ewt, roundW), computeEdgeCost_VEC(x2, y2, x4, y4, ewt, roundW));
            }

            __m256 offsetVec = _mm256_sub_ps(altEdgeWgt, solEdgeWgt); // value is negative if altEdgeWgt is better

            // compare current offset with best offsets
            __m256 mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);

            // set new bests if any
            bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
            bestIDsVec = _mm256_blendv_epi8(bestIDsVec, idsVec, _mm256_castps_si256(mask));

            // increment ids vec
            idsVec = _mm256_add_epi32(idsVec, increment);
        }

        float vecStore[AVX_VEC_SIZE];
        int idsVecStore[AVX_VEC_SIZE];

        // update the best variables
        _mm256_storeu_ps(vecStore, bestOffsetVec);
        _mm256_storeu_si256((__m256i *)idsVecStore, bestIDsVec);

        for (int i = 0; i < AVX_VEC_SIZE; i++)
        {
            if (bestFix.costOffset > vecStore[i])
            {
                bestFix.costOffset = vecStore[i];
                bestFix.edge0 = edge0;
                bestFix.edge1 = idsVecStore[i];
            }
        }
    }

    return bestFix;
}
