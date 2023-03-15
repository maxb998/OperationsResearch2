#include "tsp.h"


// MAX_THWORKSPACE MUST BE A MULTIPLE OF 4 (the 64 is chosen arbitrarly)
#define MAX_THWORSKPACE 4 * 64

typedef struct
{
    Instance *d;
    pthread_mutex_t mutex;
    size_t nextRow;
    //size_t thWorkspace;
} ThreadedInstance;


//int distMatIntegerW(Instance *d);

// method to compute edges weights with multiple threads using AVX2 instructions
static void * computeSqDistMatThread(void * d);





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
    //if (d->params.edgeWeightType > 1)
    //    LOG (LOG_LVL_ERROR, "Distance Matrix Computation: Edge weight type unsupported");

    /*  uncomment method on top also 
    if (Type requires so)
       return distMatIntegerW;
    */

    // allocate memory
    size_t n = d->nodesCount;
    //d->edgeCostMat = aligned_alloc(32, n * n * sizeof(double));
    d->edgeCost.rowSizeMem = n + AVX_VEC_SIZE - n % AVX_VEC_SIZE;   // way easier to work with avx instructions if each row start pointer is alligned
    d->edgeCost.mat = aligned_alloc(32, d->edgeCost.rowSizeMem * n * sizeof(double));

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
static void * computeSqDistMatThread(void* arg)
{
    // correctly identify data
    ThreadedInstance *th = (ThreadedInstance*)arg;

    size_t n = th->d->nodesCount;

    int pthreadMutexErrorID = 0;  // check for mutex errors
    while ( (pthreadMutexErrorID = pthread_mutex_lock(&th->mutex) == 0) && (th->nextRow < th->d->nodesCount) )    // lock mutex before checking nextRow
    {
        // CRITICAL REGION STARTED(INCLUDING WHILE CONDITION)
        // here thread gets it's workspace

        size_t row = th->nextRow;
        th->nextRow++; // prevent other threads threads to access this thread working row

        pthreadMutexErrorID = pthread_mutex_unlock (&th->mutex);
        // CRITICAL REGION END

        // check errors with mutex
        if (pthreadMutexErrorID != 0) break;

        // now the thread can compute the distance matrix inside it's workspace (row)
        register __m256d x1, y1, xDiff, yDiff, dist;

        x1 = _mm256_set1_pd(th->d->X[row]);
        y1 = _mm256_set1_pd(th->d->Y[row]);

        for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
        {
            xDiff = _mm256_sub_pd(x1, _mm256_load_pd(&th->d->X[i]));
            yDiff = _mm256_sub_pd(y1, _mm256_load_pd(&th->d->Y[i]));

            // square and sum
            dist = _mm256_add_pd(_mm256_mul_pd(xDiff, xDiff), _mm256_mul_pd(yDiff, yDiff));

            // store result in memory
            _mm256_store_pd(&th->d->edgeCost.mat[row * th->d->edgeCost.rowSizeMem + i], dist);
        }
        
        // set last elements of the row (the ones that prevents a lot of headackes with avx) to max value of double(so it never gets picked in greedy algorithms)
        for (size_t i = n; i < th->d->edgeCost.rowSizeMem; i++)
            th->d->edgeCost.mat[row * th->d->edgeCost.rowSizeMem + i] = DOUBLE_MAX;
        
    }

    // unlock mutex that was locked in the while condition (assumes that)
    pthreadMutexErrorID = pthread_mutex_unlock(&th->mutex);
    
    // Log mutex errors
    if (pthreadMutexErrorID != 0)
        LOG(LOG_LVL_CRITICAL, "computeSqDistMatThread(): Error on mutex lock/unlock with code %d. Do not trust results", pthreadMutexErrorID);

    return NULL;
}
    