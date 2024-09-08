#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

// in seconds
#define LOG_INTERVAL 30
//#define USE_MOVE_14_52_36

typedef struct
{
    float costOffset;
    int edge0;
    int edge1;
    int edge2;
} __attribute__((__aligned__(64))) _3optMoveData;

typedef struct
{
    Solution *sol;
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        float *X;
        float *Y;
    #endif
    float *costCache;
    int *sectionCopy;
    int iter;

    pthread_mutex_t mutex;
    pthread_cond_t waitUpdate;

    int nThreads;
    int threadsWaiting;
    int nextEdge;
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
    bool approxSearch;
    #endif
    bool notFinished;
    double printTimeSec;

    _3optMoveData bestFixes[MAX_THREADS];
} _3optData;

// describe 3opt move type. edge1=(1,2) edge2=(3,4) edge3=(5,6)
enum _3optMoveType {
    _3OPT_MOVE_NONE=0,
    _3OPT_MOVE_13_24=1, // 2opt move
    _3OPT_MOVE_13_25_46=2,
    _3OPT_MOVE_14_53_26=4,
    _3OPT_MOVE_15_42_36=8,
    #ifdef USE_MOVE_14_52_36
        _3OPT_MOVE_14_52_36=16
    #endif
};

typedef struct
{
    _3optData *data;
    int threadID;
} wrapper;


// Decides whether to print LOG_LVL_NOTICE benchmarking information(n° of iterations and iter/sec) at the end of the run -> Used because when a metaheuristic calls 3Opt a lot those lines really clutter a lot the console
static bool printPerformanceLog = false;
// Set global varaible
void set3OptPerformanceBenchmarkLogMT(bool val)
{
    printPerformanceLog = val;
}

static void *run3OptThread(void* arg);

// Perform solution update accordingly (invert part of the solution between selected indexes(edge0,edge1) of the bestFix)
static inline bool updateSolution(_3optData *data, _3optMoveData bestFix);

// Invert a section of the path of the solution adapting all data in data
static inline void invertSection(_3optData *data, int firstEdgePos, int lastEdgePos);

// Swap two sections inside baseAddr using backupArray to copy one of those sections
static inline void swapSections(int *baseAddr, int *backupArray, int smallSectionSize, int bigSectionSize, int smallSectionStartPos, int bigSectionStartPos, bool smallSectionIsFirst);

// Search for best possible 3Opt move in data->sol
static inline void _3OptBestFix(wrapper *w, int edge0);
#if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
    // Search for best possible 3Opt move in data->sol using approximated cost computation (needed a second function in order to have some kind of performance improvement)
    static inline void _3OptBestFixApprox(wrapper *w, int edge0);
#endif



void apply3OptBestFixMT(Solution *sol)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    sol->indexPath[n] = sol->indexPath[0];

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))

        float *costCache = malloc((n + AVX_VEC_SIZE) * 3 * sizeof(float) + n * sizeof(int));
        if (costCache == NULL)
            throwError("apply3OptBestFix: Failed to allocate memory");
        float *X = &costCache[n + AVX_VEC_SIZE];
        float *Y = &X[n + AVX_VEC_SIZE];
        int *sectionCopy = (int*)&Y[n + AVX_VEC_SIZE];

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

    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)

        float *costCache = malloc((n + AVX_VEC_SIZE) * sizeof(float) + n + sizeof(int));
        if (costCache == NULL)
            throwError("apply3OptBestFix: Failed to allocate memory");
        int *sectionCopy = (int*)&costCache[n + AVX_VEC_SIZE];

    #endif

    // build cost cache
    for (int i = 0; i < n; i++)
        #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            costCache[i] = computeEdgeCost(X[i], Y[i], X[i + 1], Y[i + 1], inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            costCache[i] = inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i+1]];
        #endif

    for (int i = n + 1; i < n + AVX_VEC_SIZE; i++) // fill remaining slots with non-interfering values
        costCache[i] = INFINITY;

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        apply3OptBestFix_fastIterativelyMT(sol, X, Y, costCache, sectionCopy);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        apply3OptBestFix_fastIterativelyMT(sol, costCache, sectionCopy);
    #endif

    free(costCache);
}

#if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
void apply3OptBestFix_fastIterativelyMT(Solution *sol, float *X, float *Y, float *costCache, int *sectionCopy)
#elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
void apply3OptBestFix_fastIterativelyMT(Solution *sol, float *costCache, int *sectionCopy)
#endif
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // check solution correspondence with X and Y when debugging
    #ifdef DEBUG
        Instance *inst = sol->instance;
        int n = inst->nNodes;
        
        if (!checkSolution(sol))
            throwError("apply3OptBestFix: Input solution is not valid");
        
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            for (int i = 0; i <= n; i++)
                if ((inst->X[sol->indexPath[i]] != X[i]) || (inst->Y[sol->indexPath[i]] != Y[i]))
                    throwError("apply3OptBestFix_fastIteratively: input mismatch between inst.X/Y[%d] = [%f, %f] and X/Y[indexPath[%d]] = [%f, %f]", i, inst->X[sol->indexPath[i]], X[i], i, inst->Y[sol->indexPath[i]], Y[i]);
        #endif

        // check costCache
        for (int i = 0; i < n; i++)
            #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
                if ((costCache[i] != computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i + 1]], inst->Y[sol->indexPath[i + 1]], inst)) && !((inst->params.mode == MODE_TABU) && (costCache[i] == -INFINITY)))
                    throwError("apply3OptBestFix_fastIteratively: input not valid: costCache isn't coherent with solution at position %d", i);
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                if ((costCache[i] != inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i+1]]) && !((inst->params.mode == MODE_TABU) && (costCache[i] == -INFINITY)))
                    throwError("apply3OptBestFix_fastIteratively: input not valid: costCache isn't coherent with solution at position %d", i);
            #endif
    #endif

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        _3optData data = { .sol=sol, .X=X, .Y=Y, .costCache=costCache, .sectionCopy=sectionCopy, .iter=0, .nThreads=sol->instance->params.nThreads, .threadsWaiting=0, .nextEdge=0, .approxSearch=true, .notFinished=true, .printTimeSec=startTime };
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        _3optData data = { .sol=sol, .X=X, .Y=Y, .costCache=costCache, .sectionCopy=sectionCopy, .iter=0, .nThreads=sol->instance->params.nThreads, .threadsWaiting=0, .nextEdge=0, .notFinished=true, .printTimeSec=startTime };
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        _3optData data = { .sol=sol, .costCache=costCache, .sectionCopy=sectionCopy, .iter=0, .nThreads=sol->instance->params.nThreads, .threadsWaiting=0, .nextEdge=0, .notFinished=true, .printTimeSec=startTime };
    #endif

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
        pthread_create(&threads[i], NULL, run3OptThread, (void*)&wrappers[i]);
    }
    
    for (int i = 0; i < data.nThreads; i++)
        pthread_join(threads[i], NULL);
    
    pthread_mutex_destroy(&data.mutex);
    pthread_cond_destroy(&data.waitUpdate);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double elapsed = cvtTimespec2Double(timeStruct) - startTime;
    if (printPerformanceLog)
    {
        LOG(LOG_LVL_NOTICE, "Total number of iterations: %lu", data.iter);
        LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)data.iter/elapsed);
    }

    sol->execTime += elapsed;
}

static void *run3OptThread(void* arg)
{
    wrapper *w = (wrapper*)arg;
    _3optData *data = w->data;

    struct timespec timeStruct;

    while (data->notFinished) // runs 3opt until no more moves are made in one iteration of 3opt
    {
        int edge0 = data->nextEdge;
        if (edge0 == data->sol->instance->nNodes)
        {
            pthread_mutex_lock(&data->mutex);
            if (edge0 == data->sol->instance->nNodes)
            {
                if (data->threadsWaiting == data->nThreads-1)
                {
                    _3optMoveData bestFix = data->bestFixes[0];
                    for (int i = 1; i < data->nThreads; i++)
                        if (data->bestFixes[i].costOffset < bestFix.costOffset)
                            bestFix = data->bestFixes[i];
                    
                    bool result = updateSolution(data, bestFix);
                    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
                        if (!result && data->approxSearch)
                        {
                            if (printPerformanceLog)
                                LOG(LOG_LVL_DEBUG, "apply3OptBestFix_fastIterativelyMT[%d]: Switching from Approximated Search to Exact Search", data->iter);
                            data->approxSearch = false;
                        }
                        else if (!result)
                            data->notFinished = false;
                    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
                        if (!result)
                            data->notFinished = false;
                    #endif

                    #ifdef DEBUG
                        if (!checkSolution(data->sol))
                            throwError("apply3OptBestFix_fastIterativelyMT: [%d] Solution is not correct", data->iter);
                    #endif
                    data->threadsWaiting = 0;
                    data->nextEdge = 0;

                    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
                    double currentTime = cvtTimespec2Double(timeStruct);
                    if (printPerformanceLog && (currentTime - data->printTimeSec > LOG_INTERVAL))
                    {   
                        LOG(LOG_LVL_INFO, "3Opt running: cost is %lf at iteration %4lu with last optimization of %lf", cvtCost2Double(data->sol->cost), data->iter, -bestFix.costOffset);
                        data->printTimeSec = currentTime;
                    }
                    data->iter++;

                    for (int i = 0; i < data->nThreads; i++)
                        data->bestFixes[i].costOffset = 0;

                    pthread_cond_broadcast(&data->waitUpdate);
                }
                else
                {
                    data->threadsWaiting++;
                    pthread_cond_wait(&data->waitUpdate, &data->mutex);
                }
            }
            else
                data->nextEdge++;
            pthread_mutex_unlock(&data->mutex);
        }
        else
            data->nextEdge++;
        if (!data->notFinished)
            break;

        #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
            if (data->approxSearch)
                _3OptBestFixApprox(w, edge0);
            else
                _3OptBestFix(w, edge0);
        #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
            _3OptBestFix(w, edge0);
        #endif
    }

    return NULL;
}

static inline bool updateSolution(_3optData *data, _3optMoveData bestFix)
{
    Solution *sol = data->sol;
    Instance *inst = sol->instance;

    int e0 = bestFix.edge0, e1 = bestFix.edge1, e2 = bestFix.edge2;

    if (bestFix.costOffset > -EPSILON)
        return false;

    // compute necessary costs
    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        float cost13 = computeEdgeCost(data->X[e0  ], data->Y[e0  ], data->X[e1  ], data->Y[e1  ], inst);
        float cost14 = computeEdgeCost(data->X[e0  ], data->Y[e0  ], data->X[e1+1], data->Y[e1+1], inst);
        float cost24 = computeEdgeCost(data->X[e0+1], data->Y[e0+1], data->X[e1+1], data->Y[e1+1], inst);
        float cost36 = computeEdgeCost(data->X[e1  ], data->Y[e1  ], data->X[e2+1], data->Y[e2+1], inst);
        float cost25 = computeEdgeCost(data->X[e0+1], data->Y[e0+1], data->X[e2  ], data->Y[e2  ], inst);
        float cost15 = computeEdgeCost(data->X[e0  ], data->Y[e0  ], data->X[e2  ], data->Y[e2  ], inst);
        float cost46 = computeEdgeCost(data->X[e1+1], data->Y[e1+1], data->X[e2+1], data->Y[e2+1], inst);
        float cost35 = computeEdgeCost(data->X[e1  ], data->Y[e1  ], data->X[e2  ], data->Y[e2  ], inst);
        float cost26 = computeEdgeCost(data->X[e0+1], data->Y[e0+1], data->X[e2+1], data->Y[e2+1], inst);
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        int *p = sol->indexPath;
        size_t n = inst->nNodes;
        float cost13 = inst->edgeCostMat[p[e0  ] * n + p[e1  ]];
        float cost14 = inst->edgeCostMat[p[e0  ] * n + p[e1+1]];
        float cost24 = inst->edgeCostMat[p[e0+1] * n + p[e1+1]];
        float cost36 = inst->edgeCostMat[p[e1  ] * n + p[e2+1]];
        float cost25 = inst->edgeCostMat[p[e0+1] * n + p[e2  ]];
        float cost15 = inst->edgeCostMat[p[e0  ] * n + p[e2  ]];
        float cost46 = inst->edgeCostMat[p[e1+1] * n + p[e2+1]];
        float cost35 = inst->edgeCostMat[p[e1  ] * n + p[e2  ]];
        float cost26 = inst->edgeCostMat[p[e0+1] * n + p[e2+1]];
    #endif

    // find out if best move type is selected
    float bestOffset = INFINITY;
    // _3OPT_MOVE_NONE
    float oldE0E1Cost = data->costCache[e0] + data->costCache[e1], oldE0E1E2Cost = oldE0E1Cost + data->costCache[e2];
    enum _3optMoveType bestMoveType = _3OPT_MOVE_NONE;

    // _3OPT_MOVE_13_24
    float offset = cost13 + cost24 - oldE0E1Cost;
    if (offset < bestOffset)
    {
        bestMoveType = _3OPT_MOVE_13_24;
        bestOffset = offset;
    }

    // _3OPT_MOVE_15_42_36
    offset = cost15 + cost36 + cost24 - oldE0E1E2Cost;
    if (offset < bestOffset)
    {
        bestMoveType = _3OPT_MOVE_15_42_36;
        bestOffset = offset;
    }

    // _3OPT_MOVE_13_25_46
    offset = cost13 + cost46 + cost25 - oldE0E1E2Cost;
    if (offset < bestOffset)
    {
        bestMoveType = _3OPT_MOVE_13_25_46;
        bestOffset = offset;
    }

    // _3OPT_MOVE_14_53_26
    offset = cost14 + cost35 + cost26 - oldE0E1E2Cost;
    if (offset < bestOffset)
    {
        bestMoveType = _3OPT_MOVE_14_53_26;
        bestOffset = offset;
    }

    // _3OPT_MOVE_14_52_36
    #ifdef USE_MOVE_14_52_36
        offset = cost14 + cost36 + cost25 - oldE0E1E2Cost;
        if (offset < bestOffset)
        {
            bestMoveType = _3OPT_MOVE_14_52_36;
            bestOffset = offset;
        }
    #endif

    #ifdef DEBUG
        if (bestOffset != bestFix.costOffset)
            LOG(LOG_LVL_WARN, "3OptBestFix -> updateSolution[%d]: recomputed cost is different from the one specified (%f vs %f)", data->iter, bestOffset, bestFix.costOffset);
    #endif

    // (check and) update cost
    float altE0=0, altE1=0, altE2=0;
    switch (bestMoveType)
    {
    case _3OPT_MOVE_NONE:
        throwError("3Opt: Error during solution update, found no matching costOffset");
        break;
    case _3OPT_MOVE_13_24:
        altE0 = cost13; altE1 = cost24; altE2 = data->costCache[e2];
        break;
    case _3OPT_MOVE_13_25_46:
        altE0 = cost13; altE1 = cost25; altE2 = cost46;
        break;
    case _3OPT_MOVE_14_53_26:
        altE0 = cost14; altE1 = cost35; altE2 = cost26;
        break;
    case _3OPT_MOVE_15_42_36:
        altE0 = cost15; altE1 = cost24; altE2 = cost36;
        break;
    #ifdef USE_MOVE_14_52_36
        case _3OPT_MOVE_14_52_36:
            altE0 = cost14; altE1 = cost25; altE2 = cost36;
            break;
    #endif
    }

    #if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
        float recompOffset = altE0 + altE1 + altE2;
        recompOffset -= (data->costCache[e0] + data->costCache[e1] + data->costCache[e2]);
        if (recompOffset > -EPSILON)
            return false;
    #endif

    sol->cost += cvtFloat2Cost(altE0) + cvtFloat2Cost(altE1) + cvtFloat2Cost(altE2) - cvtFloat2Cost(data->costCache[e0]) - cvtFloat2Cost(data->costCache[e1]) - cvtFloat2Cost(data->costCache[e2]);

    LOG(LOG_LVL_TRACE, "3Opt: [%d] Updating solution by switching edges (%d,%d,%d) with moveType %d improving cost by %f. New Cost = %lf", data->iter,
        e0, e1, e2, bestMoveType,
        bestFix.costOffset, cvtCost2Double(sol->cost));

    // invert first section if necessary
    if (bestMoveType & (_3OPT_MOVE_13_24 | _3OPT_MOVE_13_25_46 | _3OPT_MOVE_14_53_26))
        invertSection(data, e0, e1);
    // invert second section if necessary
    if (bestMoveType & (_3OPT_MOVE_13_25_46 | _3OPT_MOVE_15_42_36))
        invertSection(data, e1, e2);

    // swap sections if necessary
    #ifdef USE_MOVE_14_52_36
        if (bestFix.movetype & (_3OPT_MOVE_14_53_26 | _3OPT_MOVE_15_42_36 | _3OPT_MOVE_14_52_36))
    #else
        if (bestMoveType & (_3OPT_MOVE_14_53_26 | _3OPT_MOVE_15_42_36 ))
    #endif
    {
        // find smallest section in order to use copy that one into the sectionCopy
        int smallSectionSize = e1 - e0, smallSectionStartPos = e0 + 1;
        int bigSectionSize = e2 - e1, bigSectionStartPos = e1 + 1;
        bool smallSectionIsFirst = true;

        if (smallSectionSize > bigSectionSize)
        {
            swapElems(smallSectionSize, bigSectionSize)
            swapElems(smallSectionStartPos, bigSectionStartPos)
            smallSectionIsFirst = false;
        }

        swapSections(data->sol->indexPath, data->sectionCopy, smallSectionSize, bigSectionSize, smallSectionStartPos, bigSectionStartPos, smallSectionIsFirst);
        swapSections((int*)data->costCache, data->sectionCopy, smallSectionSize, bigSectionSize, smallSectionStartPos, bigSectionStartPos, smallSectionIsFirst);
        #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            swapSections((int*)data->X, data->sectionCopy, smallSectionSize, bigSectionSize, smallSectionStartPos, bigSectionStartPos, smallSectionIsFirst);
            swapSections((int*)data->Y, data->sectionCopy, smallSectionSize, bigSectionSize, smallSectionStartPos, bigSectionStartPos, smallSectionIsFirst);
        #endif
    }

    // update costCache with new weights
    if (bestMoveType & _3OPT_MOVE_13_24)
    {
        data->costCache[e0] = altE0;
        data->costCache[e1] = altE1;
    }
    #ifdef USE_MOVE_14_52_36
        else if (bestMoveType & (_3OPT_MOVE_14_53_26 | _3OPT_MOVE_15_42_36 | _3OPT_MOVE_14_52_36))
    #else
        else if (bestMoveType & (_3OPT_MOVE_14_53_26 | _3OPT_MOVE_15_42_36))
    #endif
    {
        data->costCache[e0] = altE0;
        data->costCache[e0 + (e2 - e1)] = altE1;
        data->costCache[e2] = altE2;
    }
    else if (bestMoveType & _3OPT_MOVE_13_25_46)
    {
        data->costCache[e0] = altE0;
        data->costCache[e1] = altE1;
        data->costCache[e2] = altE2;
    }

    #ifdef USE_MOVE_14_52_36
        if (bestFix.movetype & _3OPT_MOVE_14_52_36)
            LOG(LOG_LVL_WARN, "Finally found moveType: _3OPT_MOVE_14_52_36!");
    #endif

    return true;
}

static inline void invertSection(_3optData *data, int firstEdgePos, int lastEdgePos)
{
    int first = firstEdgePos+1, last = lastEdgePos;

    while (first < last)
    {
        swapElems(data->sol->indexPath[first], data->sol->indexPath[last])

        #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
            swapElems(data->X[first], data->X[last])
            swapElems(data->Y[first], data->Y[last])
        #endif

        first++;
        last--;
    }

    first = firstEdgePos + 1;
    last = lastEdgePos - 1;

    while (first < last)
    {
        swapElems(data->costCache[first], data->costCache[last])

        first++;
        last--;
    }
}

static inline void swapSections(int *baseAddr, int *backupArray, int smallSectionSize, int bigSectionSize, int smallSectionStartPos, int bigSectionStartPos, bool smallSectionIsFirst)
{
    for (int i = 0; i < smallSectionSize; i++)
        backupArray[i] = baseAddr[smallSectionStartPos + i];

    if (smallSectionIsFirst)
    {
        for (int i = 0; i < bigSectionSize; i++)
            baseAddr[smallSectionStartPos + i] = baseAddr[bigSectionStartPos + i];
        for (int i = 0; i < smallSectionSize; i++)
            baseAddr[smallSectionStartPos + bigSectionSize + i] = backupArray[i];
    }
    else
    {
        for (int i = bigSectionSize-1; i >= 0; i--)
            baseAddr[bigSectionStartPos + smallSectionSize + i] = baseAddr[bigSectionStartPos + i];
        for (int i = 0; i < smallSectionSize; i++)
            baseAddr[bigSectionStartPos + i] = backupArray[i];
    }
}

#if (COMPUTATION_TYPE == COMPUTE_OPTION_AVX)
static inline void _3OptBestFix(wrapper *w, int edge0)
{
    _3optData *data = w->data;
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    float *X = data->X, *Y = data->Y;

    __m256 x1 = _mm256_broadcast_ss(&X[edge0]), y1 = _mm256_broadcast_ss(&Y[edge0]);
    __m256 x2 = _mm256_broadcast_ss(&X[edge0 + 1]), y2 = _mm256_broadcast_ss(&Y[edge0 + 1]);

    for (int edge1 = 1 + edge0; edge1 < n; edge1++)
    {
        __m256 bestOffsetVec = _mm256_set1_ps(INFINITY);
        __m256i edge2Vec = _mm256_add_epi32((_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)), _mm256_set1_epi32(edge1+2));
        __m256i bestEdge2Vec = _mm256_set1_epi32(-1);


        __m256 x3 = _mm256_broadcast_ss(&X[edge1]), y3 = _mm256_broadcast_ss(&Y[edge1]);
        __m256 x4 = _mm256_broadcast_ss(&X[edge1 + 1]), y4 = _mm256_broadcast_ss(&Y[edge1 + 1]);

        float sumE0E1Cost = data->costCache[edge0] + data->costCache[edge1];

        // compute here some costs that will be needed always in order to avoid useless sqrt operations
        __m256 cost13 = computeEdgeCost_VEC(x1, y1, x3, y3, inst);
        __m256 cost14 = computeEdgeCost_VEC(x1, y1, x4, y4, inst);
        __m256 cost24 = computeEdgeCost_VEC(x2, y2, x4, y4, inst);

        // 2opt check between edge0 and edge1
        {
            __m256 offsetVec = _mm256_sub_ps(_mm256_add_ps(cost13, cost24), _mm256_set1_ps(sumE0E1Cost));
            __m256 mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
            bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
        }

        for (int edge2Pos = 2 + edge1; edge2Pos < n; edge2Pos+=AVX_VEC_SIZE)
        {
            __m256 x5 = _mm256_loadu_ps(&data->X[edge2Pos]), y5 = _mm256_loadu_ps(&data->Y[edge2Pos]);
            __m256 x6 = _mm256_loadu_ps(&data->X[edge2Pos+1]), y6 = _mm256_loadu_ps(&data->Y[edge2Pos+1]);

            // sum of the cost of the three edges in the solution
            __m256 solutionCostVec = _mm256_add_ps(_mm256_set1_ps(sumE0E1Cost), _mm256_loadu_ps(&data->costCache[edge2Pos]));

            __m256 cost36 = computeEdgeCost_VEC(x3, y3, x6, y6, inst);
            __m256 cost25 = computeEdgeCost_VEC(x2, y2, x5, y5, inst);

            //_3OPT_MOVE_14_52_36
            #ifdef USE_MOVE_14_52_36
                __m256 offsetVec = _mm256_add_ps(_mm256_add_ps(cost14, cost36), cost25);
                __m256 mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));
            #endif

            if (edge0 + 1 != edge1)
            {
                //_3OPT_MOVE_15_42_36
                __m256 offsetVec = _mm256_sub_ps(_mm256_add_ps(_mm256_add_ps(cost24, cost36), computeEdgeCost_VEC(x1, y1, x5, y5, inst)), solutionCostVec);
                __m256 mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));

                //_3OPT_MOVE_13_25_46
                offsetVec = _mm256_sub_ps(_mm256_add_ps(_mm256_add_ps(cost13, cost25), computeEdgeCost_VEC(x4, y4, x6, y6, inst)), solutionCostVec);
                mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));

                //_3OPT_MOVE_14_53_26
                offsetVec = _mm256_sub_ps(_mm256_add_ps(_mm256_add_ps(cost14, computeEdgeCost_VEC(x3, y3, x5, y5, inst)), computeEdgeCost_VEC(x2, y2, x6, y6, inst)), solutionCostVec);
                mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));
            }

            edge2Vec = _mm256_add_epi32(edge2Vec, _mm256_set1_epi32(AVX_VEC_SIZE));
        }

        float bestOffsetVecStore[AVX_VEC_SIZE];
        int bestEdge2VecStore[AVX_VEC_SIZE];

        // update the best variables
        _mm256_storeu_ps(bestOffsetVecStore, bestOffsetVec);
        _mm256_storeu_si256((__m256i *)bestEdge2VecStore, bestEdge2Vec);

        for (int i = 0; i < AVX_VEC_SIZE; i++)
        {
            if (data->bestFixes[w->threadID].costOffset > bestOffsetVecStore[i])
            {
                data->bestFixes[w->threadID].costOffset = bestOffsetVecStore[i];
                data->bestFixes[w->threadID].edge0 = edge0;
                data->bestFixes[w->threadID].edge1 = edge1;
                data->bestFixes[w->threadID].edge2 = bestEdge2VecStore[i];
            }
        }
    }
}
static inline void _3OptBestFixApprox(wrapper *w, int edge0)
{
    _3optData *data = w->data;
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    float *X = data->X, *Y = data->Y;

    __m256 x1 = _mm256_broadcast_ss(&X[edge0]), y1 = _mm256_broadcast_ss(&Y[edge0]);
    __m256 x2 = _mm256_broadcast_ss(&X[edge0 + 1]), y2 = _mm256_broadcast_ss(&Y[edge0 + 1]);

    for (int edge1 = 1 + edge0; edge1 < n; edge1++)
    {
        __m256 bestOffsetVec = _mm256_set1_ps(INFINITY);
        __m256i edge2Vec = _mm256_add_epi32((_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)), _mm256_set1_epi32(edge1+2));
        __m256i bestEdge2Vec = _mm256_set1_epi32(-1);


        __m256 x3 = _mm256_broadcast_ss(&X[edge1]), y3 = _mm256_broadcast_ss(&Y[edge1]);
        __m256 x4 = _mm256_broadcast_ss(&X[edge1 + 1]), y4 = _mm256_broadcast_ss(&Y[edge1 + 1]);

        float sumE0E1Cost = data->costCache[edge0] + data->costCache[edge1];

        // compute here some costs that will be needed always in order to avoid useless sqrt operations
        __m256 cost13 = computeApproxEdgeCost_VEC(x1, y1, x3, y3, inst);
        __m256 cost14 = computeApproxEdgeCost_VEC(x1, y1, x4, y4, inst);
        __m256 cost24 = computeApproxEdgeCost_VEC(x2, y2, x4, y4, inst);

        // 2opt check between edge0 and edge1
        {
            __m256 offsetVec = _mm256_sub_ps(_mm256_add_ps(cost13, cost24), _mm256_set1_ps(sumE0E1Cost));
            __m256 mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
            bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
        }

        for (int edge2Pos = 2 + edge1; edge2Pos < n; edge2Pos+=AVX_VEC_SIZE)
        {
            __m256 x5 = _mm256_loadu_ps(&data->X[edge2Pos]), y5 = _mm256_loadu_ps(&data->Y[edge2Pos]);
            __m256 x6 = _mm256_loadu_ps(&data->X[edge2Pos+1]), y6 = _mm256_loadu_ps(&data->Y[edge2Pos+1]);

            // sum of the cost of the three edges in the solution
            __m256 solutionCostVec = _mm256_add_ps(_mm256_set1_ps(sumE0E1Cost), _mm256_loadu_ps(&data->costCache[edge2Pos]));

            __m256 cost36 = computeApproxEdgeCost_VEC(x3, y3, x6, y6, inst);
            __m256 cost25 = computeApproxEdgeCost_VEC(x2, y2, x5, y5, inst);

            //_3OPT_MOVE_14_52_36
            #ifdef USE_MOVE_14_52_36
                __m256 offsetVec = _mm256_add_ps(_mm256_add_ps(cost14, cost36), cost25);
                __m256 mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));
            #endif

            if (edge0 + 1 != edge1)
            {
                //_3OPT_MOVE_15_42_36
                __m256 offsetVec = _mm256_sub_ps(_mm256_add_ps(_mm256_add_ps(cost24, cost36), computeApproxEdgeCost_VEC(x1, y1, x5, y5, inst)), solutionCostVec);
                __m256 mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));

                //_3OPT_MOVE_13_25_46
                offsetVec = _mm256_sub_ps(_mm256_add_ps(_mm256_add_ps(cost13, cost25), computeApproxEdgeCost_VEC(x4, y4, x6, y6, inst)), solutionCostVec);
                mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));

                //_3OPT_MOVE_14_53_26
                offsetVec = _mm256_sub_ps(_mm256_add_ps(_mm256_add_ps(cost14, computeApproxEdgeCost_VEC(x3, y3, x5, y5, inst)), computeApproxEdgeCost_VEC(x2, y2, x6, y6, inst)), solutionCostVec);
                mask = _mm256_cmp_ps(offsetVec, bestOffsetVec, _CMP_LT_OQ);
                bestOffsetVec = _mm256_blendv_ps(bestOffsetVec, offsetVec, mask);
                bestEdge2Vec = _mm256_blendv_epi8(bestEdge2Vec, edge2Vec, _mm256_castps_si256(mask));
            }

            edge2Vec = _mm256_add_epi32(edge2Vec, _mm256_set1_epi32(AVX_VEC_SIZE));
        }

        float bestOffsetVecStore[AVX_VEC_SIZE];
        int bestEdge2VecStore[AVX_VEC_SIZE];

        // update the best variables
        _mm256_storeu_ps(bestOffsetVecStore, bestOffsetVec);
        _mm256_storeu_si256((__m256i *)bestEdge2VecStore, bestEdge2Vec);

        for (int i = 0; i < AVX_VEC_SIZE; i++)
        {
            if (data->bestFixes[w->threadID].costOffset > bestOffsetVecStore[i])
            {
                data->bestFixes[w->threadID].costOffset = bestOffsetVecStore[i];
                data->bestFixes[w->threadID].edge0 = edge0;
                data->bestFixes[w->threadID].edge1 = edge1;
                data->bestFixes[w->threadID].edge2 = bestEdge2VecStore[i];
            }
        }
    }
}
#elif ((COMPUTATION_TYPE == COMPUTE_OPTION_BASE) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX))
static inline void _3OptBestFix(wrapper *w, int edge0)
{
    _3optData *data = w->data;
    Solution *sol = data->sol;
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        float *X = data->X, *Y = data->Y;
    #endif

    _3optMoveData m = {.edge0=edge0};
    
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
        float x1 = X[m.edge0], y1 = Y[m.edge0];
        float x2 = X[m.edge0 + 1], y2 = Y[m.edge0 + 1];
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        int p1 = sol->indexPath[m.edge0];
        int p2 = sol->indexPath[m.edge0 + 1];
    #endif

    for (m.edge1 = 1 + m.edge0, m.edge2=-1; m.edge1 < n; m.edge1++)
    {
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            float x3 = X[m.edge1], y3 = Y[m.edge1];
            float x4 = X[m.edge1 + 1], y4 = Y[m.edge1 + 1];
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            int p3 = sol->indexPath[m.edge1];
            int p4 = sol->indexPath[m.edge1 + 1];
        #endif

        float sumE0E1Cost = data->costCache[m.edge0] + data->costCache[m.edge1];

        // compute(fetch) here some costs that will be needed always in order to avoid useless sqrt(fetch) operations
        #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
            float cost13 = computeEdgeCost(x1, y1, x3, y3, inst);
            float cost14 = computeEdgeCost(x1, y1, x4, y4, inst);
            float cost24 = computeEdgeCost(x2, y2, x4, y4, inst);
        #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
            float cost13 = inst->edgeCostMat[p1 * (size_t)n + p3];
            float cost14 = inst->edgeCostMat[p1 * (size_t)n + p4];
            float cost24 = inst->edgeCostMat[p2 * (size_t)n + p4];
        #endif

        // 2opt check between edge0 and edge1
        m.costOffset = cost13 + cost24 - sumE0E1Cost;
        if (m.costOffset < data->bestFixes[w->threadID].costOffset)
            data->bestFixes[w->threadID] = m;

        for (m.edge2 = 2 + m.edge1; m.edge2 < n; m.edge2++)
        {
            #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                float x5 = data->X[m.edge2], y5 = data->Y[m.edge2];
                float x6 = data->X[m.edge2+1], y6 = data->Y[m.edge2+1];
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                int p5 = sol->indexPath[m.edge2];
                int p6 = sol->indexPath[m.edge2 + 1];
            #endif

            // sum of the cost of the three edges in the solution
            float costInSolution = sumE0E1Cost + data->costCache[m.edge2];

            #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                float cost36 = computeEdgeCost(x3, y3, x6, y6, inst);
                float cost25 = computeEdgeCost(x2, y2, x5, y5, inst);
            #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                float cost36 = inst->edgeCostMat[p3 * (size_t)n + p6];
                float cost25 = inst->edgeCostMat[p2 * (size_t)n + p5];
            #endif

            //_3OPT_MOVE_14_52_36
            #ifdef USE_MOVE_14_52_36
                m.costOffset = cost14 + cost36 + cost25;
                if (m.costOffset < bestFix.costOffset)
                    bestFix = m;
            #endif

            if (m.edge0 + 1 != m.edge1)
            {
                //_3OPT_MOVE_15_42_36
                #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                    m.costOffset = cost24 + cost36 + computeEdgeCost(x1, y1, x5, y5, inst) - costInSolution;
                #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                    m.costOffset = cost24 + cost36 + inst->edgeCostMat[p1 * (size_t)n + p5] - costInSolution;
                #endif
                if (m.costOffset < data->bestFixes[w->threadID].costOffset)
                    data->bestFixes[w->threadID] = m;

                //_3OPT_MOVE_13_25_46
                #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                    m.costOffset = cost13 + cost25 + computeEdgeCost(x4, y4, x6, y6, inst) - costInSolution;
                #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                    m.costOffset = cost13 + cost25 + inst->edgeCostMat[p4 * (size_t)n + p6] - costInSolution;
                #endif
                if (m.costOffset < data->bestFixes[w->threadID].costOffset)
                    data->bestFixes[w->threadID] = m;

                //_3OPT_MOVE_14_53_26
                #if (COMPUTATION_TYPE == COMPUTE_OPTION_BASE)
                    m.costOffset = cost14 + computeEdgeCost(x3, y3, x5, y5, inst) + computeEdgeCost(x2, y2, x6, y6, inst) - costInSolution;
                #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
                    m.costOffset = cost14 + inst->edgeCostMat[p3 * (size_t)n + p5] + inst->edgeCostMat[p2 * (size_t)n + p6] - costInSolution;
                #endif
                if (m.costOffset < data->bestFixes[w->threadID].costOffset)
                    data->bestFixes[w->threadID] = m;
            }
        }
    }
}
#endif
