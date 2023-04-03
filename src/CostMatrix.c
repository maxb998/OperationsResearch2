#include "CostMatrix.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

typedef struct
{
    Instance *inst;
    pthread_mutex_t mutex;
    size_t nextRow; // only value that needs critical region, used to select the row on which the thread will work on
} ThreadsData;

// method to compute edges weights with multiple threads using AVX2 instructions
static void * computeDistMatThread(void * arg);

void printCostMatrix(Instance *inst)
{
    LOG(LOG_LVL_LOG, "Printing the distance matrix");

    for (size_t row = 0; row < inst->nNodes; row++)
    {
        printf(" ");
        for (size_t col = 0; col < inst->nNodes; col++)
            printf("%.3e ", inst->edgeCostMat[row * inst->nNodes + col]);
        printf("\n");
    }
}


double computeCostMatrix(Instance *inst)
{
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    // first check that distance type is supported
    if (inst->params.edgeWeightType > ATT)
        throwError(inst, NULL, "Distance Matrix Computation: Edge weight type unsupported");

    // allocate memory
    inst->edgeCostMat = malloc(inst->nNodes * inst->nNodes * sizeof(float));
        
    // init data structure to pass to threads
    ThreadsData thInst = { .inst = inst, .nextRow = 0 };
    pthread_mutex_init (&thInst.mutex, NULL);

    // each thread compute one row of the distance matrix before calling mutex to get a new workspace

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < inst->params.nThreads; i++)
    {
        pthread_create (&threads[i], NULL, computeDistMatThread, &thInst);
        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld CREATED", i);
    }

    for (size_t i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join (threads[i], NULL);
        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld FINISHED", i);
    }

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    double elapsed = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);
    return elapsed;
}

// EUCLIDEAN DISTANCE
static void * computeDistMatThread(void* arg)
{
    // correctly identify data
    ThreadsData *th = (ThreadsData*)arg;
    Instance *inst = th->inst;

    size_t n = inst->nNodes;

    // initialize mask outside loop to avoid doing it over and over(useless beacuse the compiler should do this automatically but whatever)
    int mask[AVX_VEC_SIZE];

    // since it has the same value for each row we can define it here, outside the loops
    size_t lastVecElemsCount = (n % AVX_VEC_SIZE);
    if (lastVecElemsCount == 0)
        lastVecElemsCount = 8;
    for (size_t j = 0; j < lastVecElemsCount; j++)
        mask[j] = -1; // this row last elements
    for (size_t j = lastVecElemsCount; j < AVX_VEC_SIZE; j++)
        mask[j] = 0; // next row first elements (DO NOT WRITE THIS IN THE STORE -> so we use maskstore)

    while ( (pthread_mutex_lock(&th->mutex) == 0) && (th->nextRow < th->inst->nNodes) )    // lock mutex before checking nextRow
    {
        // CRITICAL REGION STARTED(INCLUDING WHILE CONDITION) ######################
        // here thread gets it's workspace

        size_t row = th->nextRow;
        th->nextRow++; // prevent other threads threads to access this thread working row

        pthread_mutex_unlock (&th->mutex);
        // CRITICAL REGION END #####################################################

        // now the thread can compute the distance matrix inside it's workspace (row)
        __m256 x1 = _mm256_set1_ps(th->inst->X[row]);
        __m256 y1 = _mm256_set1_ps(th->inst->Y[row]);

        size_t i;
        for (i = 0; i < n - AVX_VEC_SIZE; i += AVX_VEC_SIZE)
        {
            __m256 x2 = _mm256_loadu_ps(&th->inst->X[i]);
            __m256 y2 = _mm256_loadu_ps(&th->inst->Y[i]);

            __m256 dist = computeEdgeCost_VEC(x1, y1, x2, y2, inst->params.edgeWeightType, inst->params.roundWeights);

            // store result in memory
            _mm256_storeu_ps(&inst->edgeCostMat[row * n + i], dist);
        }
        
        // last elements must be done differently since a vector may override the contents of the next line

        // do the same as in the loop
        __m256 x2 = _mm256_loadu_ps(&th->inst->X[i]);
        __m256 y2 = _mm256_loadu_ps(&th->inst->Y[i]);
        __m256 dist = computeEdgeCost_VEC(x1, y1, x2, y2, inst->params.edgeWeightType, inst->params.roundWeights);

        // maskstore the result avoiding to overwrite data
        _mm256_maskstore_ps(&inst->edgeCostMat[row * n + i], _mm256_loadu_si256((__m256i*)mask), dist);
    }

    // unlock mutex that was locked in the while condition
    pthread_mutex_unlock(&th->mutex);

    pthread_exit(NULL);
}
    