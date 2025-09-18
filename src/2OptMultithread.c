#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

//#define DEBUG

// in seconds
#define LOG_INTERVAL 30

typedef struct
{
    float costOffset;
    int edge0;
    int edge1;
} __attribute__((__aligned__(64))) _2optMoveData;

typedef struct
{
    Solution *sol;
    float *X;
    float *Y;
    float *costCache;
    int iter;

    pthread_mutex_t mutex;
    pthread_cond_t waitUpdate;

    int nThreads;
    int threadsWaiting;
    int nextEdge;
    bool approxSearch;
    bool notFinished;
    double printTimeSec;

    _2optMoveData bestFixes[MAX_THREADS];
} _2optData;

typedef struct
{
    _2optData *data;
    int threadID;
} wrapper;


// Decides whether to print LOG_LVL_NOTICE benchmarking information(n° of iterations and iter/sec) at the end of the run -> Used because when a metaheuristic calls 2opt a lot those lines really clutter a lot the console
static bool print2OptLog = false;
// Set global varaible
void set2OptPerformanceBenchmarkLogMT(bool val)
{
    print2OptLog = val;
}


// Perform solution update accordingly (invert part of the solution between selected indexes(edge0,edge1) of the bestFix)
static inline bool updateSolution(_2optData *data, _2optMoveData bestFix);

// Search for best possible 2opt move in data->sol
static inline void _2OptBestFixBase(wrapper *w, int edge0);
static inline void _2OptBestFixMatrix(wrapper *w, int edge0);
static inline void _2OptBestFixAVX(wrapper *w, int edge0);

static inline void _2OptBestFixApproxAVX(wrapper *w, int edge0);

static void *run2OptThread(void* arg);


void apply2OptBestFixMT(Solution *sol)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    float *costCache = NULL, *X = NULL, *Y = NULL;

    sol->indexPath[n] = sol->indexPath[0];

    if (inst->params.compType & (COMP_BASE|COMP_AVX))
    {
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

        // build cost cache
        for (int i = 0; i < n; i++)
            costCache[i] = computeEdgeCost(X[i], Y[i], X[i + 1], Y[i + 1], inst);
    }
    else
    {
        costCache = malloc((n + AVX_VEC_SIZE) * sizeof(float));
        if (costCache == NULL)
            throwError("apply2OptBestFix: Failed to allocate memory");
        
        // build cost cache
        for (int i = 0; i < n; i++)
            costCache[i] = inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i+1]];
    }

    for (int i = n + 1; i < n + AVX_VEC_SIZE; i++) // fill remaining slots with non-interfering values
        costCache[i] = INFINITY;

    apply2OptBestFix_fastIterativelyMT(sol, X, Y, costCache);

    free(costCache);
}

int apply2OptBestFix_fastIterativelyMT(Solution *sol, float *X, float *Y, float *costCache)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    Instance *inst = sol->instance;
    
    // check solution correspondence with X and Y when debugging
    #ifdef DEBUG
        int n = inst->nNodes;
        
        if (!checkSolution(sol))
            throwError("apply2OptBestFix: Input solution is not valid");
        
        if (inst->params.compType & COMP_AVX)
        {
            for (int i = 0; i <= n; i++)
                if ((inst->X[sol->indexPath[i]] != X[i]) || (inst->Y[sol->indexPath[i]] != Y[i]))
                    throwError("apply2OptBestFix_fastIterativelyMT: input mismatch between inst.X/Y[%d] = [%f, %f] and X/Y[indexPath[%d]] = [%f, %f]", i, inst->X[sol->indexPath[i]], X[i], i, inst->Y[sol->indexPath[i]], Y[i]);
        }

        // check costCache
        if (inst->params.compType & (COMP_BASE|COMP_AVX))
        {
            for (int i = 0; i < n; i++)
                if ((costCache[i] != computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i + 1]], inst->Y[sol->indexPath[i + 1]], inst)) && !((inst->params.mode == MODE_TABU) && (costCache[i] == -INFINITY)))
                    throwError("apply2OptBestFix_fastIterativelyMT: input not valid: costCache isn't coherent with solution at position %d", i);
        }
        else
        {
            for (int i = 0; i < n; i++)
                if ((costCache[i] != inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i+1]]) && !((inst->params.mode == MODE_TABU) && (costCache[i] == -INFINITY)))
                    throwError("apply2OptBestFix_fastIterativelyMT: input not valid: costCache isn't coherent with solution at position %d", i);
        }
    #endif
        
    _2optData data = { .sol=sol, .X=X, .Y=Y, .costCache=costCache, .iter=0, .nThreads=sol->instance->params.nThreads, .threadsWaiting=0, .nextEdge=0, .approxSearch=false, .notFinished=true, .printTimeSec=startTime };
    if (inst->params.compType & COMP_AVX)
        data.approxSearch = true;

    for (int i = 0; i < data.nThreads; i++)
        data.bestFixes[i].costOffset = 0;

    pthread_mutex_init(&data.mutex, NULL);
    pthread_cond_init(&data.waitUpdate, NULL);
    
    pthread_t threads[MAX_THREADS];
    wrapper wrappers[MAX_THREADS];
    for (int i = 0; i < data.nThreads; i++)
    {
        wrappers[i].data = &data;
        wrappers[i].threadID = i;
        pthread_create(&threads[i], NULL, run2OptThread, (void*)&wrappers[i]);
    }
    
    for (int i = 0; i < data.nThreads; i++)
        pthread_join(threads[i], NULL);
    
    pthread_mutex_destroy(&data.mutex);
    pthread_cond_destroy(&data.waitUpdate);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double elapsed = cvtTimespec2Double(timeStruct) - startTime;
    if (print2OptLog)
    {
        LOG(LOG_LVL_NOTICE, "Total number of iterations: %lu", data.iter);
        LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)data.iter/elapsed);
    }

    sol->execTime += elapsed;
    return data.iter-1;
}

static void *run2OptThread(void* arg)
{
    wrapper *w = (wrapper*)arg;
    _2optData *data = w->data;
    Instance *inst = data->sol->instance;

    struct timespec timeStruct;

    while (data->notFinished) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        int edge0 = data->nextEdge;
        if (edge0 >= data->sol->instance->nNodes)
        {
            pthread_mutex_lock(&data->mutex);
            if (data->threadsWaiting == data->nThreads-1)
            {
                _2optMoveData bestFix = data->bestFixes[0];
                for (int i = 1; i < data->nThreads; i++)
                    if (data->bestFixes[i].costOffset < bestFix.costOffset)
                        bestFix = data->bestFixes[i];
                
                bool result = updateSolution(data, bestFix);
                if (inst->params.compType & COMP_AVX)
                {
                    if (!result && data->approxSearch)
                    {
                        if (print2OptLog)
                            LOG(LOG_LVL_DEBUG, "apply2OptBestFix_fastIterativelyMT[%d]: Switching from Approximated Search to Exact Search", data->iter);
                        data->approxSearch = false;
                    }
                    else if (!result)
                        data->notFinished = false;
                }
                else
                {
                    if (!result)
                        data->notFinished = false;
                }

                #ifdef DEBUG
                    if (!checkSolution(data->sol))
                        throwError("apply2OptBestFix_fastIterativelyMT: [%d] Solution is not correct", data->iter);
                #endif
                data->threadsWaiting = 0;
                data->nextEdge = 0;

                clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
                double currentTime = cvtTimespec2Double(timeStruct);
                if (print2OptLog && (currentTime - data->printTimeSec > LOG_INTERVAL))
                {   
                    LOG(LOG_LVL_INFO, "2Opt running: cost is %lf at iteration %4lu with last optimization of %lf", cvtCost2Double(data->sol->cost), data->iter, -bestFix.costOffset);
                    data->printTimeSec = currentTime;
                }
                data->iter++;

                pthread_cond_broadcast(&data->waitUpdate);
            }
            else
            {
                data->threadsWaiting++;
                pthread_cond_wait(&data->waitUpdate, &data->mutex);
            }
            data->bestFixes[w->threadID].costOffset = 0;
            pthread_mutex_unlock(&data->mutex);
        }
        else
            data->nextEdge++;
        if (!data->notFinished)
            break;
        
        switch (inst->params.compType)
        {
        case COMP_BASE:
            _2OptBestFixBase(w, edge0);
            break;
        case COMP_MATRIX:
            _2OptBestFixMatrix(w, edge0);
            break;
        case COMP_AVX:
            if (data->approxSearch)
                _2OptBestFixApproxAVX(w, edge0);
            else
                _2OptBestFixAVX(w, edge0);
            break;
        default:
            break;
        }
    }

    return NULL;
}

static inline bool updateSolution(_2optData *data, _2optMoveData bestFix)
{
    Solution *sol = data->sol;
    Instance *inst = sol->instance;

    //float oldEdge0Cost, oldEdge1Cost;
    float altEdge0Cost, altEdge1Cost;

    if (inst->params.compType & (COMP_BASE|COMP_AVX))
    {
        altEdge0Cost = computeEdgeCost(data->X[bestFix.edge0], data->Y[bestFix.edge0], data->X[bestFix.edge1], data->Y[bestFix.edge1], inst);
        altEdge1Cost = computeEdgeCost(data->X[bestFix.edge0+1], data->Y[bestFix.edge0+1], data->X[bestFix.edge1+1], data->Y[bestFix.edge1+1], inst);
    }
    else
    {
        int *indexPath = sol->indexPath;
        int n = inst->nNodes;
        altEdge0Cost = inst->edgeCostMat[(size_t)indexPath[bestFix.edge0] * (size_t)n + (size_t)indexPath[bestFix.edge1]];
        altEdge1Cost = inst->edgeCostMat[(size_t)indexPath[bestFix.edge0+1] * (size_t)n + (size_t)indexPath[bestFix.edge1+1]];
    }

    bestFix.costOffset = altEdge0Cost + altEdge1Cost;
    bestFix.costOffset -= (data->costCache[bestFix.edge0] + data->costCache[bestFix.edge1]);
    if ((bestFix.costOffset > -EPSILON) || (bestFix.edge0 == -1))
        return false;

    // update cost
    sol->cost += cvtFloat2Cost(altEdge0Cost) + cvtFloat2Cost(altEdge1Cost) - cvtFloat2Cost(data->costCache[bestFix.edge0]) - cvtFloat2Cost(data->costCache[bestFix.edge1]);

    if (print2OptLog)
    {
        LOG(LOG_LVL_TRACE, "2Opt: [%d] Updating solution by switching edge (%d,%d)[%d] with edge (%d,%d)[%d] reducing cost by %f. New Cost = %lf", data->iter,
            sol->indexPath[bestFix.edge0], sol->indexPath[bestFix.edge0 + 1], bestFix.edge0,
            sol->indexPath[bestFix.edge1], sol->indexPath[bestFix.edge1 + 1], bestFix.edge1,
            -bestFix.costOffset, cvtCost2Double(sol->cost));
    }

    for (int s = bestFix.edge0+1, b = bestFix.edge1; s < b; s++, b--)
        swapElems(sol->indexPath[s], sol->indexPath[b])
    
    if (inst->params.compType & (COMP_BASE|COMP_AVX))
    {
        for (int s = bestFix.edge0+1, b = bestFix.edge1; s < b; s++, b--)
        {
            swapElems(data->X[s], data->X[b])
            swapElems(data->Y[s], data->Y[b])
        }
    }

    // update cost cache
    data->costCache[bestFix.edge0] = altEdge0Cost;
    data->costCache[bestFix.edge1] = altEdge1Cost;

    for (int s = bestFix.edge0+1, b = bestFix.edge1-1; s < b; s++, b--)
        swapElems(data->costCache[s], data->costCache[b])

    return true;
}

static inline void _2OptBestFixBase(wrapper *w, int edge0)
{
    _2optData *data = w->data;
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    _2optMoveData currFix = { .costOffset=0, .edge0=edge0 };

    float partSolEdgeWgt = data->costCache[currFix.edge0];

    for (currFix.edge1 = 2 + currFix.edge0; (currFix.edge1 < n - 1) || ((currFix.edge1 < n) && (currFix.edge0 > 0)); currFix.edge1++)
    {
        float solEdgeWgt = partSolEdgeWgt + data->costCache[currFix.edge1];

        // check the combined weight other combination of edges
        float altEdgeWgt = computeEdgeCost(data->X[currFix.edge0], data->Y[currFix.edge0], data->X[currFix.edge1], data->Y[currFix.edge1], inst) + 
                            computeEdgeCost(data->X[currFix.edge0 + 1], data->Y[currFix.edge0 + 1], data->X[currFix.edge1 + 1], data->Y[currFix.edge1 + 1], inst);

        currFix.costOffset = altEdgeWgt - solEdgeWgt;
        // update local best if current one is better
        if (data->bestFixes[w->threadID].costOffset > currFix.costOffset)
            data->bestFixes[w->threadID] = currFix;
    }
}

static inline void _2OptBestFixMatrix(wrapper *w, int edge0)
{
    _2optData *data = w->data;
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    _2optMoveData currFix = { .costOffset=0, .edge0=edge0 };

    float partSolEdgeWgt = data->costCache[currFix.edge0];

    for (currFix.edge1 = 2 + currFix.edge0; (currFix.edge1 < n - 1) || ((currFix.edge1 < n) && (currFix.edge0 > 0)); currFix.edge1++)
    {
        float solEdgeWgt = partSolEdgeWgt + data->costCache[currFix.edge1];

        // check the combined weight other combination of edges
        float altEdgeWgt = inst->edgeCostMat[(size_t)sol->indexPath[currFix.edge0] * (size_t)n + (size_t)sol->indexPath[currFix.edge1]] + 
                            inst->edgeCostMat[(size_t)sol->indexPath[currFix.edge0 + 1] * (size_t)n + (size_t)sol->indexPath[currFix.edge1 + 1]];

        currFix.costOffset = altEdgeWgt - solEdgeWgt;
        // update local best if current one is better
        if (data->bestFixes[w->threadID].costOffset > currFix.costOffset)
            data->bestFixes[w->threadID] = currFix;
    }
}


static inline void _2OptBestFixAVX(wrapper *w, int edge0)
{
    _2optData *data = w->data;
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    float *X = data->X, *Y = data->Y;

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
        if (data->bestFixes[w->threadID].costOffset > vecStore[i])
        {
            data->bestFixes[w->threadID].costOffset = vecStore[i];
            data->bestFixes[w->threadID].edge0 = edge0;
            data->bestFixes[w->threadID].edge1 = idsVecStore[i];
        }
    }
}

static inline void _2OptBestFixApproxAVX(wrapper *w, int edge0)
{
    _2optData *data = w->data;
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    float *X = data->X, *Y = data->Y;

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
            altEdgeWgt = _mm256_add_ps(computeApproxEdgeCost_VEC(x1, y1, x3, y3, inst), computeApproxEdgeCost_VEC(x2, y2, x4, y4, inst));
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
        if (data->bestFixes[w->threadID].costOffset > vecStore[i])
        {
            data->bestFixes[w->threadID].costOffset = vecStore[i];
            data->bestFixes[w->threadID].edge0 = edge0;
            data->bestFixes[w->threadID].edge1 = idsVecStore[i];
        }
    }
}


