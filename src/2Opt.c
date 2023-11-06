#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

//#define DEBUG

// in seconds
#define LOG_INTERVAL 30

typedef struct
{
    Solution *sol;
    float *X;
    float *Y;
    float *costCache;
    int iter;
} _2optData;


typedef struct
{
    float costOffset;
    int edge0;
    int edge1;
} _2optMoveData;


// Decides whether to print LOG_LVL_NOTICE benchmarking information(n° of iterations and iter/sec) at the end of the run -> Used because when a metaheuristic calls 2opt a lot those lines really clutter a lot the console
static bool printPerformanceLog = false;
// Set global varaible
void set2OptPerformanceBenchmarkLog(bool val)
{
    printPerformanceLog = val;
}


// Perform solution update accordingly (invert part of the solution between selected indexes(edge0,edge1) of the bestFix)
static inline void updateSolution(_2optData *data, _2optMoveData bestFix);

// Search for best possible 2opt move in data->sol
static inline _2optMoveData _2OptBestFix(_2optData *data);


void apply2OptBestFix(Solution *sol)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    float *X = NULL, *Y = NULL, *costCache = NULL;

    sol->indexPath[n] = sol->indexPath[0];

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)

        costCache = malloc((n + AVX_VEC_SIZE) * 3 * sizeof(float));
        if (costCache == NULL)
            throwError("apply2OptBestFix: Failed to allocate memory");
        X = &costCache[n + AVX_VEC_SIZE];
        Y = &X[n + AVX_VEC_SIZE];

        for (int i = 0; i < n; i++) // fill X and Y
        {
            X[i] = inst->X[sol->indexPath[i]];
            Y[i] = inst->Y[sol->indexPath[i]];
        }
        X[n] = X[0];
        Y[n] = Y[0];

        for (int i = n + 1; i < n + AVX_VEC_SIZE; i++) // fill remaining slots with non-interfering values
        {
            X[i] = INFINITY;
            Y[i] = INFINITY;
        }

    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))

        costCache = malloc((n + AVX_VEC_SIZE) * sizeof(float));
        if (costCache == NULL)
            throwError("apply2OptBestFix: Failed to allocate memory");

    #endif

    // build cost cache
    for (int i = 0; i < n; i++)
        #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            costCache[i] = computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i + 1]], inst->Y[sol->indexPath[i + 1]], inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            costCache[i] = inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i+1]];
        #endif

    for (int i = n + 1; i < n + AVX_VEC_SIZE; i++) // fill remaining slots with non-interfering values
        costCache[i] = INFINITY;

    apply2OptBestFix_fastIteratively(sol, X, Y, costCache);

    free(costCache);
}

void apply2OptBestFix_fastIteratively(Solution *sol, float *X, float *Y, float *costCache)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);
    double printTimeSec = startTime;

    // check solution correspondence with X and Y when debugging
    #ifdef DEBUG
        Instance *inst = sol->instance;
        int n = inst->nNodes;
        
        if (!checkSolution(sol))
            throwError("apply2OptBestFix: Input solution is not valid");
        
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            for (int i = 0; i <= n; i++)
                if ((inst->X[sol->indexPath[i]] != X[i]) || (inst->Y[sol->indexPath[i]] != Y[i]))
                    throwError("apply2OptBestFix_fastIteratively: input mismatch between inst.X/Y[%d] = [%f, %f] and X/Y[indexPath[%d]] = [%f, %f]", i, inst->X[sol->indexPath[i]], X[i], i, inst->Y[sol->indexPath[i]], Y[i]);
        #endif

        // check costCache
        for (int i = 0; i < n; i++)
            #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                if (costCache[i] != computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i + 1]], inst->Y[sol->indexPath[i + 1]], inst->params.edgeWeightType, inst->params.roundWeights))
                    throwError("apply2OptBestFix_fastIteratively: input not valid: costCache isn't coherent with solution at position %d", i);
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                if (costCache[i] != inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i+1]])
                    throwError("apply2OptBestFix_fastIteratively: input not valid: costCache isn't coherent with solution at position %d", i);
            #endif
    #endif

    _2optData data = { .sol=sol, .X=X, .Y=Y, .costCache=costCache, .iter=0 };

    bool notFinishedFlag = true;
    while (notFinishedFlag) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        _2optMoveData bestFix = _2OptBestFix(&data);

        if (bestFix.costOffset < -EPSILON)
            updateSolution(&data, bestFix);
        else
            notFinishedFlag = false;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        double currentTime = cvtTimespec2Double(timeStruct);
        if (printPerformanceLog && (currentTime - printTimeSec > LOG_INTERVAL))
        {   
            LOG(LOG_LVL_LOG, "2Opt running: cost is %lf at iteration %4lu with last optimization of %lf", cvtCost2Double(sol->cost), data.iter, -bestFix.costOffset);
            printTimeSec = currentTime;
        }
        data.iter++;
    }

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double elapsed = cvtTimespec2Double(timeStruct) - startTime;
    if (printPerformanceLog)
    {
        LOG(LOG_LVL_NOTICE, "Total number of iterations: %lu", data.iter);
        LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)data.iter/elapsed);
    }

    sol->execTime += elapsed;
}

static inline void updateSolution(_2optData *data, _2optMoveData bestFix)
{
    Solution *sol = data->sol;
    Instance *inst = sol->instance;

    //float oldEdge0Cost, oldEdge1Cost;
    float altEdge0Cost, altEdge1Cost;

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        altEdge0Cost = computeEdgeCost(data->X[bestFix.edge0], data->Y[bestFix.edge0], data->X[bestFix.edge1], data->Y[bestFix.edge1], inst);
        altEdge1Cost = computeEdgeCost(data->X[bestFix.edge0+1], data->Y[bestFix.edge0+1], data->X[bestFix.edge1+1], data->Y[bestFix.edge1+1], inst);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        int *indexPath = sol->indexPath;
        altEdge0Cost = computeEdgeCost(inst->X[indexPath[bestFix.edge0]], inst->Y[indexPath[bestFix.edge0]], inst->X[indexPath[bestFix.edge1]], inst->Y[indexPath[bestFix.edge1]], inst);
        altEdge1Cost = computeEdgeCost(inst->X[indexPath[bestFix.edge0+1]], inst->Y[indexPath[bestFix.edge0+1]], inst->X[indexPath[bestFix.edge1+1]], inst->Y[indexPath[bestFix.edge1+1]], inst);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        int *indexPath = sol->indexPath;
        int n = inst->nNodes;
        altEdge0Cost = inst->edgeCostMat[(size_t)indexPath[bestFix.edge0] * (size_t)n + (size_t)indexPath[bestFix.edge1]];
        altEdge1Cost = inst->edgeCostMat[(size_t)indexPath[bestFix.edge0+1] * (size_t)n + (size_t)indexPath[bestFix.edge1+1]];
    #endif

    // update cost
    sol->cost += cvtFloat2Cost(altEdge0Cost) + cvtFloat2Cost(altEdge1Cost) - cvtFloat2Cost(data->costCache[bestFix.edge0]) - cvtFloat2Cost(data->costCache[bestFix.edge1]);

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

    LOG(LOG_LVL_EVERYTHING, "2Opt: [%d] Updating solution by switching edge (%d,%d) with edge (%d,%d) improving cost by %f. New Cost = %lf", data->iter,
        sol->indexPath[bestFix.edge0], sol->indexPath[bestFix.edge0 + 1],
        sol->indexPath[bestFix.edge1], sol->indexPath[bestFix.edge1 + 1],
        bestFix.costOffset, cvtCost2Double(sol->cost));

    while (smallID < bigID)
    {
        swapElems(sol->indexPath[smallID], sol->indexPath[bigID])
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            swapElems(data->X[smallID], data->X[bigID])
            swapElems(data->Y[smallID], data->Y[bigID])
        #endif

        smallID++;
        bigID--;
    }

    // update cost cache
    data->costCache[bestFix.edge0] = altEdge0Cost;
    data->costCache[bestFix.edge1] = altEdge1Cost;

    smallID = bestFix.edge0 + 1;
    bigID = bestFix.edge1 - 1;

    while (smallID < bigID)
    {
        swapElems(data->costCache[smallID], data->costCache[bigID])

        smallID++;
        bigID--;
    }
}

#if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
static inline _2optMoveData _2OptBestFix(_2optData *data)
{
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    float *X = data->X, *Y = data->Y;

    _2optMoveData bestFix = { .costOffset=0 };

    for (int edge0 = 0; edge0 < n - 1; edge0++) // check for one edge at a time every other edge(except already checked)
    {

        __m256 x1 = _mm256_broadcast_ss(&X[edge0]), y1 = _mm256_broadcast_ss(&Y[edge0]);
        __m256 x2 = _mm256_broadcast_ss(&X[edge0 + 1]), y2 = _mm256_broadcast_ss(&Y[edge0 + 1]);
        __m256 partialSolEdgeWgt = _mm256_broadcast_ss(&data->costCache[edge0]);
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
                solEdgeWgt = _mm256_add_ps(partialSolEdgeWgt, _mm256_loadu_ps(&data->costCache[i]));

                altEdgeWgt = _mm256_add_ps(computeEdgeCost_VEC(x1, y1, x3, y3, inst), computeEdgeCost_VEC(x2, y2, x4, y4, inst));
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
#elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
static inline _2optMoveData _2OptBestFix(_2optData *data)
{
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    _2optMoveData bestFix = { .costOffset=0 };

    for (_2optMoveData currFix = { .edge0=0 }; currFix.edge0 < n - 1; currFix.edge0++) // check for one edge at a time every other edge(except already checked)
    {
        float partSolEdgeWgt = data->costCache[currFix.edge0];

        for (currFix.edge1 = 2 + currFix.edge0; (currFix.edge1 < n - 1) || ((currFix.edge1 < n) && (currFix.edge0 > 0)); currFix.edge1++)
        {
            float solEdgeWgt = partSolEdgeWgt + data->costCache[currFix.edge1];

            // check the combined weight other combination of edges
            float altEdgeWgt;
            #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                altEdgeWgt = computeEdgeCost(inst->X[sol->indexPath[currFix.edge0]], inst->Y[sol->indexPath[currFix.edge0]], inst->X[sol->indexPath[currFix.edge1]], inst->Y[sol->indexPath[currFix.edge1]], inst) + 
                             computeEdgeCost(inst->X[sol->indexPath[currFix.edge0 + 1]], inst->Y[sol->indexPath[currFix.edge0 + 1]], inst->X[sol->indexPath[currFix.edge1 + 1]], inst->Y[sol->indexPath[currFix.edge1 + 1]], inst);
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                altEdgeWgt = inst->edgeCostMat[(size_t)sol->indexPath[currFix.edge0] * (size_t)n + (size_t)sol->indexPath[currFix.edge1]] + 
                             inst->edgeCostMat[(size_t)sol->indexPath[currFix.edge1 + 1] * (size_t)n + (size_t)sol->indexPath[currFix.edge0 + 1]];
            #endif

            currFix.costOffset = altEdgeWgt - solEdgeWgt;
            // update local best if current one is better
            if (bestFix.costOffset > currFix.costOffset)
                bestFix = currFix;
        }
    }

    return bestFix;
}
#endif

