#include "tsp.h"

typedef struct
{
    Instance *d;
    pthread_mutex_t mutex;
    size_t nextRow; // only value that needs critical region, used to select the row on which the thread will work on
    char fatal;
} ThreadedInstance;

// method to compute edges weights with multiple threads using AVX2 instructions
static void * computeDistMatThread(void * d);


// #############################################################################################################
// EXACT DISTANCE FUNCTIONS DEFINITIONS

static inline __m256 euclideanCostSquared2D(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 xDiff, yDiff, dist;

    xDiff = _mm256_sub_ps(x1, x2);
    yDiff = _mm256_sub_ps(y1, y2);
    dist = _mm256_add_ps(_mm256_mul_ps(xDiff, xDiff), _mm256_mul_ps(yDiff, yDiff));

    return dist;
}

// https://stackoverflow.com/questions/63599391/find-absolute-in-avx
static inline __m256 manhattanCost2D(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 xDiff, yDiff, dist, absMask = _mm256_set1_ps(-0.0F);
    
    xDiff = _mm256_sub_ps(x1, x2);
    yDiff = _mm256_sub_ps(y1, y2);
    xDiff = _mm256_andnot_ps(absMask, xDiff);
    yDiff = _mm256_andnot_ps(absMask, yDiff);
    dist = _mm256_add_ps(xDiff, yDiff);

    return dist;
}

static inline __m256 maximumCost2D(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 xDiff, yDiff, dist;

    xDiff = _mm256_sub_ps(x1, x2);
    yDiff = _mm256_sub_ps(y1, y2);
    dist = _mm256_max_ps(xDiff, yDiff);

    return dist;
}

static inline __m256 attCostSquared2D(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 xDiff, yDiff, dist, vec10 = _mm256_set1_ps(10.F);

    xDiff = _mm256_sub_ps(x1, x2);
    yDiff = _mm256_sub_ps(y1, y2);
    dist = _mm256_add_ps(_mm256_mul_ps(xDiff, xDiff), _mm256_mul_ps(yDiff, yDiff));
    dist = _mm256_div_ps(dist, vec10);

    return dist;
}

// ###########################################################################################################
// ROUNDED DISTANCE FUNCTIONS DEFINITIONS

static inline __m256i euclideanCostSquared2DRounded(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    return _mm256_cvttps_epi32(_mm256_sqrt_ps(euclideanCostSquared2D(x1,y1,x2,y2)));
}

static inline __m256i manhattanCost2DRounded(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    return _mm256_cvttps_epi32(manhattanCost2D(x1,y1,x2,y2));
}

static inline __m256i maximumCost2DRounded(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    return _mm256_cvttps_epi32(maximumCost2D(x1,y1,x2,y2));
}

static inline __m256i attCostSquared2DRounded(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 roundedDist = attCostSquared2D(x1,y1,x2,y2);
    roundedDist = _mm256_ceil_ps(_mm256_sqrt_ps(roundedDist));
    return _mm256_cvttps_epi32(roundedDist);
}
            
// #############################################################################################################






void printDistanceMatrix(Instance *d, int showEndRowPlaceholder)
{
    size_t rowSizeToPrint;
    if (showEndRowPlaceholder == 1)
        rowSizeToPrint = d->edgeCost.rowSizeMem;
    else
        rowSizeToPrint = d->nodesCount;

    LOG(LOG_LVL_LOG, "Printing the distance matrix");

    if (d->params.roundWeights == 0)
    {
        for (size_t row = 0; row < d->nodesCount; row++)
        {
            printf(" ");
            for (size_t col = 0; col < rowSizeToPrint; col++)
                printf("%.2e ", d->edgeCost.mat[row * d->edgeCost.rowSizeMem + col]);
            printf("\n");
        }
    }
    else
    {
        for (size_t row = 0; row < d->nodesCount; row++)
        {
            printf(" ");
            for (size_t col = 0; col < rowSizeToPrint; col++)
                printf("%10d ", d->edgeCost.roundedMat[row * d->edgeCost.rowSizeMem + col]);
            printf("\n");
        }
    }
}



double computeDistanceMatrix(Instance *d)
{
    double totCpuTime = 0.;
    clock_t start, end;

    start = clock();

    // first check that distance type is supported
    if (d->params.edgeWeightType > ATT)
        throwError(d, "Distance Matrix Computation: Edge weight type unsupported");

    // allocate memory
    size_t n = d->nodesCount;

    if (n % AVX_VEC_SIZE > 0)
        d->edgeCost.rowSizeMem = n + AVX_VEC_SIZE - n % AVX_VEC_SIZE;   // way easier to work with avx instructions if each row start pointer is alligned
    else
        d->edgeCost.rowSizeMem = n;

    if (d->params.roundWeights == 0)
        d->edgeCost.mat = aligned_alloc(32, d->edgeCost.rowSizeMem * n * sizeof(float));
    else
        d->edgeCost.roundedMat = aligned_alloc(32, d->edgeCost.rowSizeMem * n * sizeof(int));
        
    // init data structure to pass to threads
    ThreadedInstance thInst = { .d = d, .nextRow = 0, .fatal = 0};
    pthread_mutex_init (&thInst.mutex, NULL);

    // each thread compute one row of the distance matrix before calling mutex to get a new workspace

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        pthread_create (&threads[i], NULL, computeDistMatThread, &thInst);
        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld CREATED", i);
    }
    
    end = clock();
    totCpuTime += ((double)(end-start)) / CLOCKS_PER_SEC;

    // wait threads to finish and get time of each thread
    double threadsTime[MAX_THREADS], *returnPtr;
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        pthread_join (threads[i], (void**)&returnPtr);
        threadsTime[i] = *returnPtr;
        free(returnPtr);

        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld finished in %.3e seconds", i, threadsTime[i]);
    }

    double threadMax = 0.;
    for (size_t i = 0; i < d->params.threadsCount; i++)
        if (threadMax < threadsTime[i])
            threadMax = threadsTime[i];

    if (thInst.fatal == 1)
        throwError(d, "Distance Matrix Computation: An error has occured while the threads were working. Check the code");

    return totCpuTime + threadMax;
}

// EUCLIDEAN DISTANCE
static void * computeDistMatThread(void* d)
{
    // measure first tick
    clock_t start, end;
    start = clock();

    // correctly identify data
    ThreadedInstance *th = (ThreadedInstance*)d;

    size_t n = th->d->nodesCount;

    int pthreadMutexErrorID = 0;  // check for mutex errors
    while ( (pthreadMutexErrorID = pthread_mutex_lock(&th->mutex) == 0) && (th->nextRow < th->d->nodesCount) )    // lock mutex before checking nextRow
    {
        // CRITICAL REGION STARTED(INCLUDING WHILE CONDITION) ######################
        // here thread gets it's workspace

        size_t row = th->nextRow;
        th->nextRow++; // prevent other threads threads to access this thread working row

        pthreadMutexErrorID = pthread_mutex_unlock (&th->mutex);
        // CRITICAL REGION END #####################################################

        // check errors with mutex
        if (pthreadMutexErrorID != 0) break;

        // now the thread can compute the distance matrix inside it's workspace (row)
        register __m256 x1, y1, x2, y2, dist;

        x1 = _mm256_set1_ps(th->d->X[row]);
        y1 = _mm256_set1_ps(th->d->Y[row]);

        if (th->d->params.roundWeights == 0)
        {
            for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
            {
                x2 = _mm256_load_ps(&th->d->X[i]);
                y2 = _mm256_load_ps(&th->d->Y[i]);

                switch (th->d->params.edgeWeightType)
                {
                case EUC_2D:
                    dist = euclideanCostSquared2D(x1, y1, x2, y2); break;

                case MAN_2D:
                    dist = manhattanCost2D(x1, y1, x2, y2); break;

                case MAX_2D:
                    dist = maximumCost2D(x1, y1, x2, y2); break;

                case ATT:
                    dist = attCostSquared2D(x1, y1, x2, y2); break;

                default:
                    LOG(LOG_LVL_CRITICAL, "Thread is trying to compute distance with unsupported edge weight type. Check the code");
                    th->fatal = 1;
                    return NULL;
                }

                // store result in memory
                _mm256_store_ps(&th->d->edgeCost.mat[row * th->d->edgeCost.rowSizeMem + i], dist);
            }
            // set last elements of the row (the ones that prevents a lot of headackes with avx) to infinity(so it never gets picked in greedy algorithms)
            for (size_t i = n; i < th->d->edgeCost.rowSizeMem; i++)
                th->d->edgeCost.mat[row * th->d->edgeCost.rowSizeMem + i] = INFINITY;
        }
        else
        {
            register __m256i roundedDist;

            for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
            {
                switch (th->d->params.edgeWeightType)
                {
                case EUC_2D:
                    roundedDist = euclideanCostSquared2DRounded(x1, y1, x2, y2); break;

                case MAN_2D:
                    roundedDist = manhattanCost2DRounded(x1, y1, x2, y2); break;

                case MAX_2D:
                    roundedDist = maximumCost2DRounded(x1, y1, x2, y2); break;

                case ATT:
                    roundedDist = attCostSquared2DRounded(x1, y1, x2, y2); break;

                default:
                    LOG(LOG_LVL_CRITICAL, "Thread is trying to compute distance with unsupported edge weight type. Check the code");
                    th->fatal = 1;
                    return NULL;
                }

                // store result in memory
                _mm256_store_ps((float*)&th->d->edgeCost.roundedMat[row * th->d->edgeCost.rowSizeMem + i], (__m256)roundedDist);
            }
            // set last elements of the row (the ones that prevents a lot of headackes with avx) to infinity(so it never gets picked in greedy algorithms)
            for (size_t i = n; i < th->d->edgeCost.rowSizeMem; i++)
                th->d->edgeCost.roundedMat[row * th->d->edgeCost.rowSizeMem + i] = INT_MAX;
        }
    }

    // unlock mutex that was locked in the while condition
    pthreadMutexErrorID = pthread_mutex_unlock(&th->mutex);
    
    // Log mutex errors
    if (pthreadMutexErrorID != 0)
    {
        LOG(LOG_LVL_ERROR, "computeSqDistMatThread(): Error on mutex lock/unlock with code %d. Do not trust results", pthreadMutexErrorID);
        th->fatal = 1;
        return NULL;
    }


    end =  clock();

    double *cpuTimeUsed = malloc(sizeof(double));
    *cpuTimeUsed = ((double)(end-start)) / CLOCKS_PER_SEC;

    pthread_exit(cpuTimeUsed);
}
    