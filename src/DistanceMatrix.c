#include "tsp.h"

typedef struct
{
    Instance *d;
    pthread_mutex_t mutex;
    size_t nextRow; // only value that needs critical region, used to select the row on which the thread will work on
} ThreadedInstance;




//int distMatIntegerW(Instance *d);


// method to compute edges weights with multiple threads using AVX2 instructions
static void * computeSqDistMatThread(void * d);






// #############################################################################################################
// DISTANCE FUNCTIONS

#define euclideanCostSquared2D(x1,y1,dist,i,X,Y) \
            static __m256 xDiff, yDiff;\
            xDiff = _mm256_sub_ps(x1, _mm256_load_ps(&X[i]));\
            yDiff = _mm256_sub_ps(y1, _mm256_load_ps(&Y[i]));\
            dist = _mm256_add_ps(_mm256_mul_ps(xDiff, xDiff), _mm256_mul_ps(yDiff, yDiff));

// https://stackoverflow.com/questions/63599391/find-absolute-in-avx
#define manhattanCost2D(x1,y1,dist,i,X,Y) \
            register __m256 xDiff, yDiff, absMask = _mm256_set1_ps(-0.0F);\
            xDiff = _mm256_sub_ps(x1, _mm256_load_ps(&X[i]));\
            yDiff = _mm256_sub_ps(y1, _mm256_load_ps(&Y[i]));\
            xDiff = _mm256_andnot_ps(xDiff, absMask);\
            yDiff = _mm256_andnot_ps(yDiff, absMask);\
            dist = _mm256_add_ps(xDiff, yDiff);

#define maximumCost2D(x1,y1,dist,i,X,Y) \
            static __m256 xDiff, yDiff; \
            xDiff = _mm256_sub_ps(x1, _mm256_load_ps(&X[i]));\
            yDiff = _mm256_sub_ps(y1, _mm256_load_ps(&Y[i]));\
            dist = _mm256_max_ps(xDiff, yDiff);

#define attCostSquared2D(x1,y1,dist,i,X,Y) \
            register __m256 xDiff, yDiff, roundedDiff, vec10 = _mm256_set1_ps(10.F); \
            xDiff = _mm256_sub_ps(x1, _mm256_load_ps(&X[i]));\
            yDiff = _mm256_sub_ps(y1, _mm256_load_ps(&Y[i]));\
            dist = _mm256_add_ps(_mm256_mul_ps(xDiff, xDiff), _mm256_mul_ps(yDiff, yDiff));\
            dist = _mm256_div_ps(dist, vec10);

// ROUNDED DISTANCES

#define euclideanCostSquared2DRounded(x1,y1,dist,roundedDist,i,X,Y) \
            euclideanCostSquared2D(x1,y1,dist,i,X,Y)\
            roundedDist = _mm256_cvttps_epi32(dist);

#define manhattanCost2DRounded(x1,y1,dist,roundedDist,i,X,Y) \
            manhattanCost2D(x1,y1,dist,i,X,Y)\
            roundedDist = _mm256_cvttps_epi32(dist);

#define maximumCost2DRounded(x1,y1,dist,roundedDist,i,X,Y) \
            maximumCost2D(x1,y1,dist,i,X,Y)\
            roundedDist = _mm256_cvttps_epi32(dist);

#define attCostSquared2DRounded(x1,y1,dist,roundedDist,i,X,Y) \
            attCostSquared2D(x1,y1,dist,i,X,Y)\
            dist = _mm256_sqrt(dist);\
            dist = _mm256_ceil_ps(dist);\
            roundedDist = _mm256_cvtps_epi32(dist);
            
// #############################################################################################################






void printDistanceMatrix(Instance *d, int showEndRowPlaceholder)
{
    size_t rowSizeToPrint;
    if (showEndRowPlaceholder == 1)
        rowSizeToPrint = d->edgeCost.rowSizeMem;
    else
        rowSizeToPrint = d->nodesCount;

    LOG(LOG_LVL_LOG, "Printing the distance matrix");

    for (size_t row = 0; row < d->nodesCount; row++)
    {
        printf(" ");
        for (size_t col = 0; col < rowSizeToPrint; col++)
            printf("%.2e ", d->edgeCost.mat[row * d->edgeCost.rowSizeMem + col]);
        printf("\n");
    }
}



int computeSquaredDistanceMatrix(Instance *d)
{
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
    ThreadedInstance thInst = { .d = d, .nextRow = 0 };
    pthread_mutex_init (&thInst.mutex, NULL);

    // each thread compute one row of the distance matrix before calling mutex to get a new workspace

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        pthread_create (&threads[i], NULL, computeSqDistMatThread, &thInst);
        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld CREATED", i);
    }
    
    // wait threads to finish
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        pthread_join (threads[i], NULL);
        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld FINISHED", i);
    }

    return 0;
}

// EUCLIDEAN DISTANCE
static void * computeSqDistMatThread(void* d)
{
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
        register __m256 x1, y1, dist;

        x1 = _mm256_set1_ps(th->d->X[row]);
        y1 = _mm256_set1_ps(th->d->Y[row]);

        for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
        {
            euclideanCostSquared2D(x1,y1,dist,i,th->d->X,th->d->Y)

            // store result in memory
            _mm256_store_ps(&th->d->edgeCost.mat[row * th->d->edgeCost.rowSizeMem + i], dist);
        }
        
        // set last elements of the row (the ones that prevents a lot of headackes with avx) to infinity(so it never gets picked in greedy algorithms)
        for (size_t i = n; i < th->d->edgeCost.rowSizeMem; i++)
            th->d->edgeCost.mat[row * th->d->edgeCost.rowSizeMem + i] = INFINITY;
        
    }

    // unlock mutex that was locked in the while condition (assumes that)
    pthreadMutexErrorID = pthread_mutex_unlock(&th->mutex);
    
    // Log mutex errors
    if (pthreadMutexErrorID != 0)
        LOG(LOG_LVL_CRITICAL, "computeSqDistMatThread(): Error on mutex lock/unlock with code %d. Do not trust results", pthreadMutexErrorID);

    return NULL;
}
    