#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

#define USE_FAST_SOLUTION_UPDATE 1
#define LOG_INTERVAL 5

static inline void updateSolutionNN(Solution *sol, int bestOffsetEdges[2], float bestOffset);

// Loads inverts and stores two vectors of floats according to update solution procedure. Speeds Up update procedure
static inline void invertVectorSizeElemsF(float *firstPtr, float *secondPtr);

// Loads inverts and stores two vectors of int according to update solution procedure. Speeds Up update procedure
static inline void invertVectorSizeElemsD(int *firstPtr, int *secondPtr);

// Basic approach that computes the costs every time, one at a time
static inline unsigned long _2optBestFixBase(Solution *sol);

// Approach that requires to have the precomputed costs
static inline unsigned long _2optBestFixCostMatrix(Solution *sol);

// Fastest approach since the mutli-threaded variant spends too much time in waiting other threads(probably wakeup from sleep costs more than one full iteration on one thread)
static inline unsigned long _2OptBestFixAVX(Solution *sol);

// Faster type of approach that computes costs every time using avx vectorization
static inline unsigned long _2OptBestFixAVXMultiThread(Solution *sol);

double apply2OptBestFix(Solution *sol, enum _2OptOptions option)
{
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    Instance *inst = sol->instance;

    // check if option is valid
    if (option < 0 || option > sizeof(enum _2OptOptions))
    {
        LOG(LOG_LVL_WARNING, "Option given to 2Opt (%d) is not valid. Using _2OPT_AVX_ST as default", option);
        option = _2OPT_AVX_ST;
    }

    // check integrity if debugging
    if (inst->params.logLevel >= LOG_LVL_DEBUG)
    {
        // always check solution
        if (!checkSolution(sol))
            throwError(inst, sol, "apply2OptBestFix: Input solution is not valid");

        if (option == _2OPT_PRECOMPUTED_COSTS)
        {
            if (inst->edgeCostMat == NULL)
            {
                LOG(LOG_LVL_WARNING, "2Opt: Required cost matrix, as specified in option parameter, is not available/not been computed.\n Switching to _2OPT_AVX");
                option = _2OPT_AVX_ST;
            }
        }
    }

    unsigned long iterNum = 0;

    switch (option)
    {
    case _2OPT_AVX_ST:
        iterNum = _2OptBestFixAVX(sol);
        break;

    case _2OPT_BASE:
        iterNum = _2optBestFixBase(sol);
        break;

    case _2OPT_PRECOMPUTED_COSTS:
        iterNum = _2optBestFixCostMatrix(sol);
        break;

    case _2OPT_AVX_MT:
        iterNum = _2OptBestFixAVXMultiThread(sol);
        break;
    }

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    double elapsed = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);
    
    LOG(LOG_LVL_NOTICE, "Total number of iterations: %lu", iterNum);
    LOG(LOG_LVL_NOTICE, "Iterations-per-second: %lf", (double)iterNum/elapsed);

    return elapsed;
}

static inline void updateSolutionNN(Solution *sol, int bestOffsetEdges[2], float bestOffset)
{
    // update cost
    sol->cost += (double)bestOffset;

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

    int smallID = bestOffsetEdges[0] + 1, bigID = bestOffsetEdges[1];

    static int updateCount = 0;
    updateCount++;

    LOG(LOG_LVL_EVERYTHING, "2Opt: [%d] Updating solution by switching edge (%d,%d) with edge (%d,%d) improving cost by %f", updateCount,
        sol->indexPath[bestOffsetEdges[0]], sol->indexPath[bestOffsetEdges[0] + 1],
        sol->indexPath[bestOffsetEdges[1]], sol->indexPath[bestOffsetEdges[1] + 1], bestOffset);

    if (USE_FAST_SOLUTION_UPDATE)
    {
        while (bigID - smallID >= 2 * AVX_VEC_SIZE) // condition checks that we two 256 bits vectors fit in between smallID and bigID and can be used to copy the data
        {
            bigID -= AVX_VEC_SIZE;

            invertVectorSizeElemsF(&sol->X[smallID], &sol->X[bigID + 1]);
            invertVectorSizeElemsF(&sol->Y[smallID], &sol->Y[bigID + 1]);

            invertVectorSizeElemsD(&sol->indexPath[smallID], &sol->indexPath[bigID + 1]);

            smallID += AVX_VEC_SIZE;
        }
    }

    while (smallID < bigID)
    {
        { // swap index path elements
            register int temp;
            swapElems(sol->indexPath[smallID], sol->indexPath[bigID], temp);
        }

        { // swap solution coordinates
            register float temp;
            swapElems(sol->X[smallID], sol->X[bigID], temp);
            swapElems(sol->Y[smallID], sol->Y[bigID], temp);
        }

        smallID++;
        bigID--;
    }
}

static inline void invertVectorSizeElemsF(float *firstPtr, float *secondPtr)
{
    static int inversionMask[AVX_VEC_SIZE] = {7, 6, 5, 4, 3, 2, 1, 0};

    // load data into vectors
    __m256 firstVec = _mm256_loadu_ps(firstPtr);
    __m256 secondVec = _mm256_loadu_ps(secondPtr);

    // invert smallVec and bigVec
    firstVec = _mm256_permutevar8x32_ps(firstVec, _mm256_loadu_si256((__m256i *)inversionMask));
    secondVec = _mm256_permutevar8x32_ps(secondVec, _mm256_loadu_si256((__m256i *)inversionMask));

    // store inverted
    _mm256_storeu_ps(firstPtr, secondVec);
    _mm256_storeu_ps(secondPtr, firstVec);
}

static inline void invertVectorSizeElemsD(int *firstPtr, int *secondPtr)
{
    static int inversionMask[AVX_VEC_SIZE] = {7, 6, 5, 4, 3, 2, 1, 0};

    // load data into vectors
    __m256i firstVec = _mm256_loadu_si256((__m256i *)firstPtr);
    __m256i secondVec = _mm256_loadu_si256((__m256i *)secondPtr);

    // invert smallVec and bigVec
    firstVec = _mm256_permutevar8x32_epi32(firstVec, _mm256_loadu_si256((__m256i *)inversionMask));
    secondVec = _mm256_permutevar8x32_epi32(secondVec, _mm256_loadu_si256((__m256i *)inversionMask));

    // store inverted
    _mm256_storeu_si256((__m256i *)firstPtr, secondVec);
    _mm256_storeu_si256((__m256i *)secondPtr, firstVec);
}

static inline unsigned long _2optBestFixBase(Solution *sol)
{
    Instance *inst = sol->instance;

    // setup "shortcuts" variables to declutter the code
    int n = inst->nNodes;
    enum EdgeWeightType edgeWgtType = inst->params.edgeWeightType;
    bool roundFlag = inst->params.roundWeights;

    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double printTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    unsigned long iterNum = 0;

    int finishedFlag = 0;
    while (finishedFlag == 0) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        // setup local values to avoid calling mutex too often
        float bestOffset = 0;
        int bestOffsetEdges[2];

        // do while there are edges to be checked
        for (int edgeID = 0; edgeID < n - 1; edgeID++) // check for one edge at a time every other edge(except already checked)
        {

            float solEdgeWgt1 = computeEdgeCost(sol->X[edgeID], sol->Y[edgeID], sol->X[edgeID + 1], sol->Y[edgeID + 1], edgeWgtType, roundFlag); /// inst->edgeCostMat[sol->indexPath[edgeID] * n + sol->indexPath[edgeID + 1]];
            for (int i = 2 + edgeID; (i < n - 1) || ((i < n) && (edgeID > 0)); i++)
            {
                float solEdgeWgt2 = computeEdgeCost(sol->X[i], sol->Y[i], sol->X[i + 1], sol->Y[i + 1], edgeWgtType, roundFlag); // inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i + 1]];

                // check the combined weight other combination of edges
                float altEdge1 = computeEdgeCost(sol->X[edgeID], sol->Y[edgeID], sol->X[i], sol->Y[i], edgeWgtType, roundFlag);                 // inst->edgeCostMat[sol->indexPath[edgeID] * n + sol->indexPath[i]];
                float altEdge2 = computeEdgeCost(sol->X[edgeID + 1], sol->Y[edgeID + 1], sol->X[i + 1], sol->Y[i + 1], edgeWgtType, roundFlag); // inst->edgeCostMat[sol->indexPath[i + 1] * n + sol->indexPath[edgeID + 1]];
                if (solEdgeWgt1 + solEdgeWgt2 > altEdge1 + altEdge2)
                {
                    // update local best if current one is better
                    float currentOffset = (altEdge1 + altEdge2) - (solEdgeWgt1 + solEdgeWgt2);
                    if (bestOffset > currentOffset)
                    {
                        bestOffset = currentOffset;
                        bestOffsetEdges[0] = edgeID;
                        bestOffsetEdges[1] = i;
                    }
                }
            }
        }

        // update solution if good enough optimization is found
        if (bestOffset < -EPSILON)
            updateSolutionNN(sol, bestOffsetEdges, bestOffset);
        else
            finishedFlag = 1;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        double currTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
        if (currTimeSec - printTimeSec > LOG_INTERVAL)
        {
            LOG(LOG_LVL_LOG, "Solution optimization in progress: cost is %lf at iteration %4lu with last optimization of %f", sol->cost, iterNum, -bestOffset);
            printTimeSec = currTimeSec;
        }
        iterNum++;
    }

    return iterNum;
}

static inline unsigned long _2optBestFixCostMatrix(Solution *sol)
{
    Instance *inst = sol->instance;

    // setup "shortcuts" variables to declutter the code
    int n = inst->nNodes;

    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double printTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    unsigned long iterNum = 0;

    int finishedFlag = 0;
    while (finishedFlag == 0) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        // setup local values to avoid calling mutex too often
        float bestOffset = 0;
        int bestOffsetEdges[2] = {0, 0};

        // do while there are edges to be checked
        for (int edgeID = 0; edgeID < n - 1; edgeID++) // check for one edge at a time every other edge(except already checked)
        {
            float solEdgeWgt1 = inst->edgeCostMat[sol->indexPath[edgeID] * n + sol->indexPath[edgeID + 1]];
            for (int i = 2 + edgeID; (i < n - 1) || ((i < n) && (edgeID > 0)); i++)
            {
                float solEdgeWgt2 = inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i + 1]];

                // check the combined weight other combination of edges
                float altEdge1 = inst->edgeCostMat[sol->indexPath[edgeID] * n + sol->indexPath[i]];
                float altEdge2 = inst->edgeCostMat[sol->indexPath[i + 1] * n + sol->indexPath[edgeID + 1]];
                if (solEdgeWgt1 + solEdgeWgt2 > altEdge1 + altEdge2)
                {
                    // update local best if current one is better
                    float currentOffset = (altEdge1 + altEdge2) - (solEdgeWgt1 + solEdgeWgt2);
                    if (bestOffset > currentOffset)
                    {
                        bestOffset = currentOffset;
                        bestOffsetEdges[0] = edgeID;
                        bestOffsetEdges[1] = i;
                    }
                }
            }
        }
        if (bestOffset < -EPSILON)
            updateSolutionNN(sol, bestOffsetEdges, bestOffset);
        else
            finishedFlag = 1;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        double currTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
        if (currTimeSec - printTimeSec > LOG_INTERVAL)
        {
            LOG(LOG_LVL_LOG, "Solution optimization in progress: cost is %lf at iteration %4lu with last optimization of %f", sol->cost, iterNum, -bestOffset);
            printTimeSec = currTimeSec;
        }
        iterNum++;
    }

    return iterNum;
}

static inline unsigned long _2OptBestFixAVX(Solution *sol)
{
    Instance *inst = sol->instance;

    // setup "shortcuts" variables to declutter the code
    int n = inst->nNodes;
    enum EdgeWeightType edgeWgtType = inst->params.edgeWeightType;
    bool roundFlag = inst->params.roundWeights;

    float vecStore[AVX_VEC_SIZE];
    int idsVecStore[AVX_VEC_SIZE];

    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double printTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    unsigned long iterNum = 0;

    int finishedFlag = 0;
    while (finishedFlag == 0) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        float bestOffset = 0;
        int bestOffsetEdges[2] = {0};

        // do while there are edges to be checked
        for (int edgeID = 0; edgeID < n - 1; edgeID++) // check for one edge at a time every other edge(except already checked)
        {

            __m256 x1 = _mm256_broadcast_ss(&sol->X[edgeID]), y1 = _mm256_broadcast_ss(&sol->Y[edgeID]);
            __m256 x2 = _mm256_broadcast_ss(&sol->X[edgeID + 1]), y2 = _mm256_broadcast_ss(&sol->Y[edgeID + 1]);
            __m256 partialSolEdgeWgt = computeEdgeCost_VEC(x1, y1, x2, y2, edgeWgtType, roundFlag);
            __m256 bestOffsetVec = _mm256_set1_ps(INFINITY);

            __m256i idsVec = _mm256_add_epi32((_mm256_set_epi32(9, 8, 7, 6, 5, 4, 3, 2)), _mm256_set1_epi32(edgeID));
            __m256i increment = _mm256_set1_epi32(AVX_VEC_SIZE);
            __m256i bestIDsVec = _mm256_set1_epi32(-1);

            for (int i = 2 + edgeID; (i < n - 1) || ((i < n) && (edgeID > 0)); i += AVX_VEC_SIZE)
            {
                __m256 altEdgeWgt, solEdgeWgt;
                { // scope "force" a thing compiler should do automatically -> x3,y3,x4,y4 destroyed as soon as we don't need them anymore
                    __m256 x3 = _mm256_loadu_ps(&sol->X[i]), y3 = _mm256_loadu_ps(&sol->Y[i]);
                    __m256 x4 = _mm256_loadu_ps(&sol->X[i + 1]), y4 = _mm256_loadu_ps(&sol->Y[i + 1]);
                    solEdgeWgt = _mm256_add_ps(partialSolEdgeWgt, computeEdgeCost_VEC(x3, y3, x4, y4, edgeWgtType, roundFlag));

                    altEdgeWgt = _mm256_add_ps(computeEdgeCost_VEC(x1, y1, x3, y3, edgeWgtType, roundFlag), computeEdgeCost_VEC(x2, y2, x4, y4, edgeWgtType, roundFlag));
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

            // update the best variables
            _mm256_storeu_ps(vecStore, bestOffsetVec);
            _mm256_storeu_si256((__m256i *)idsVecStore, bestIDsVec);

            for (int i = 0; i < AVX_VEC_SIZE; i++)
            {
                if (bestOffset > vecStore[i])
                {
                    bestOffset = vecStore[i];
                    bestOffsetEdges[0] = edgeID;
                    bestOffsetEdges[1] = idsVecStore[i];
                }
            }
        }

        if (bestOffset < -EPSILON)
            updateSolutionNN(sol, bestOffsetEdges, bestOffset);
        else
            finishedFlag = 1;

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        double currTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
        if (currTimeSec - printTimeSec >= LOG_INTERVAL) //(iterNum % (n / 50) == 0)
        {
            LOG(LOG_LVL_LOG, "Solution optimization in progress: cost is %lf at iteration %4lu with last optimization of %f", sol->cost, iterNum, -bestOffset);
            printTimeSec = currTimeSec;
        }

        iterNum++;
    }

    return iterNum;
}

// ###############################################################################################################################################################
// MULTITHREADED APPROACH HERE
// ###############################################################################################################################################################

typedef struct
{
    Solution *sol;

    // id of the next edge to be checked against all solution edges
    int nextEdge;

    // flag that is set to 1 when the bestSolution has stopped improving
    char finishedFlag;

    // best optimization cost offset stored here.
    // When negative a cost-reducing optimization has been found otherwise it's positive(or a very small negative to prevent numerical issues)
    float bestOffset;

    // best optimization configuration: ids of the nodes that give the best costOffset when "switched" with 2opt
    int bestOffsetEdges[2];

    // Each thread set threadFinishFlag[threadID] = 1 when it's waiting for other thread to finish and update the bestSolution
    volatile int threadFinishFlag[MAX_THREADS];

    // mutex
    pthread_mutex_t mutex;

    double printTimeSec;
    unsigned long iterNum;
} ThreadsData;

// struct that allows to give a meaningful and unique id to each thread
typedef struct
{
    ThreadsData *th;
    int threadID;
} ThreadsDataWithID;

// Function that each threads executes
static void *_2OptBestFixAVXThread(void *arg);

static inline unsigned long _2OptBestFixAVXMultiThread(Solution *sol)
{
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);

    ThreadsData th = {
        .sol = sol,
        .nextEdge = 0,
        .finishedFlag = 0,
        .bestOffset = 0,
        .printTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0,
        .iterNum = 0};
    Instance *inst = sol->instance;

    // init mutexes & condition
    pthread_mutex_init(&th.mutex, NULL);

    // struct necessary to give each thread it's identifier correctly
    ThreadsDataWithID thIDs[MAX_THREADS];
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        thIDs[i].th = &th;
        thIDs[i].threadID = i;
    }

    // start threads
    pthread_t workers[MAX_THREADS];
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_create(&workers[i], NULL, _2OptBestFixAVXThread, &thIDs[i]);
        LOG(LOG_LVL_DEBUG, "2-Opt: Thread %ld created", i);
    }

    // wait for threads
    for (int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join(workers[i], NULL);
        LOG(LOG_LVL_DEBUG, "2-Opt: Thread %ld finished", i);
    }

    // destroy mutexes and conditions
    pthread_mutex_destroy(&th.mutex);

    return th.iterNum;
}

static void *_2OptBestFixAVXThread(void *arg)
{
    ThreadsDataWithID *thID = (ThreadsDataWithID *)arg;
    ThreadsData *th = thID->th;
    int threadID = thID->threadID;
    Solution *sol = th->sol;
    Instance *inst = sol->instance;

    // setup "shortcuts" variables to declutter the code
    int n = inst->nNodes;
    enum EdgeWeightType edgeWgtType = inst->params.edgeWeightType;
    bool roundFlag = inst->params.roundWeights;

    float vecStore[AVX_VEC_SIZE];
    int idsVecStore[AVX_VEC_SIZE];

    while (th->finishedFlag == 0) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        // setup local values to avoid calling mutex too often
        float localBestOffset = 0;
        int localBestOffsetID[2];

        // do while there are edges to be checked
        while (pthread_mutex_trylock(&th->mutex) != 0)
            ;
        while (th->nextEdge < n - 1) // check for one edge at a time every other edge(except already checked)
        {
            int edgeID = th->nextEdge; // goes from 0 (means edge from first node in solution to second) to n which is the number of nodes
            th->nextEdge++;
            pthread_mutex_unlock(&th->mutex);

            __m256 x1 = _mm256_broadcast_ss(&sol->X[edgeID]), y1 = _mm256_broadcast_ss(&sol->Y[edgeID]);
            __m256 x2 = _mm256_broadcast_ss(&sol->X[edgeID + 1]), y2 = _mm256_broadcast_ss(&sol->Y[edgeID + 1]);
            __m256 partialSolEdgeWgt = computeEdgeCost_VEC(x1, y1, x2, y2, edgeWgtType, roundFlag); // inst->edgeCostMat[sol->indexPath[edgeID] * n + sol->indexPath[edgeID + 1]];
            __m256 bestOffsetVec = _mm256_set1_ps(INFINITY);

            __m256i idsVec = _mm256_add_epi32((_mm256_set_epi32(9, 8, 7, 6, 5, 4, 3, 2)), _mm256_set1_epi32(edgeID));
            __m256i increment = _mm256_set1_epi32(AVX_VEC_SIZE);
            __m256i bestIDsVec = _mm256_set1_epi32(-1);

            for (int i = 2 + edgeID; (i < n - 1) || ((i < n) && (edgeID > 0)); i += AVX_VEC_SIZE)
            {
                __m256 altEdgeWgt, solEdgeWgt;
                { // scope "force" a thing compiler should do automatically -> x3,y3,x4,y4 destroyed as soon as we don't need them anymore
                    __m256 x3 = _mm256_loadu_ps(&sol->X[i]), y3 = _mm256_loadu_ps(&sol->Y[i]);
                    __m256 x4 = _mm256_loadu_ps(&sol->X[i + 1]), y4 = _mm256_loadu_ps(&sol->Y[i + 1]);
                    solEdgeWgt = _mm256_add_ps(partialSolEdgeWgt, computeEdgeCost_VEC(x3, y3, x4, y4, edgeWgtType, roundFlag));

                    altEdgeWgt = _mm256_add_ps(computeEdgeCost_VEC(x1, y1, x3, y3, edgeWgtType, roundFlag), computeEdgeCost_VEC(x2, y2, x4, y4, edgeWgtType, roundFlag));
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

            // update the localBest variables
            _mm256_storeu_ps(vecStore, bestOffsetVec);
            _mm256_storeu_si256((__m256i *)idsVecStore, bestIDsVec);

            for (int i = 0; i < AVX_VEC_SIZE; i++)
            {
                if (localBestOffset > vecStore[i])
                {
                    localBestOffset = vecStore[i];
                    localBestOffsetID[0] = edgeID;
                    localBestOffsetID[1] = idsVecStore[i];
                }
            }

            while (pthread_mutex_trylock(&th->mutex) != 0)
                ;
            // pthread_mutex_lock(&th->mutex);
        }
        // When this point is reached the mutex is still locked

        // When all possibilities for one edge have been checked, update thread shared bestoffset if localbest is better
        if (th->bestOffset > localBestOffset)
        {
            th->bestOffset = localBestOffset;
            th->bestOffsetEdges[0] = localBestOffsetID[0];
            th->bestOffsetEdges[1] = localBestOffsetID[1];
        }

        th->nextEdge = 0;
        th->threadFinishFlag[threadID] = 1;

        int lastThreadFlag = 1; // if it remains one than this is the last thread going through this part at the end of the 2opt move search
        for (int i = 0; i < inst->params.nThreads; i++)
            if (th->threadFinishFlag[i] == 0)
                lastThreadFlag = 0;

        // signal updater thread to update solution only when all worker threads have finished
        if (lastThreadFlag == 1)
        {
            if (th->bestOffset < -EPSILON)
                updateSolutionNN(sol, th->bestOffsetEdges, th->bestOffset);
            else
                th->finishedFlag = 1;

            struct timespec currT;
            clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
            double currTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
            if (currTimeSec - th->printTimeSec >= LOG_INTERVAL) //(iterNum % (n / 50) == 0)
            {
                LOG(LOG_LVL_LOG, "Solution optimization in progress: cost is %lf at iteration %4lu with last optimization of %f", sol->cost, th->iterNum, -th->bestOffset);
                th->printTimeSec = currTimeSec;
            }
            th->iterNum++;

            // reset offset
            th->bestOffset = 0;

            // reset all the threadFinishFlag to 0
            for (int i = 0; i < inst->params.nThreads; i++)
                th->threadFinishFlag[i] = 0;
        }
        pthread_mutex_unlock(&th->mutex);

        while (th->threadFinishFlag[threadID] != 0);
    }

    pthread_exit(NULL);
}
