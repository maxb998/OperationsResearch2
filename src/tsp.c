#include "tsp.h"


#define MAX_THWORSKPACE 1024

typedef struct
{
    Instance *d;
    pthread_mutex_t mutex;
    size_t freeElement;
    size_t thWorkspace;
} ThreadedInstance;


//int distMatIntegerW(Instance *d);

// method to compute edges weights with multiple threads using AVX2 instructions
static void * computeSqDistMatThread(void * d);

int computeSquaredDistanceMatrix(Instance *d)
{
    // first check that distance type is supported
    if (d->params.edgeWeightType > 1)
        LOG (LOG_LVL_ERROR, "Distance Matrix Computation: Edge weight type unsupported");

    /*  uncomment method on top also 
    if (Type requires so)
       return distMatIntegerW;
    */

    // allocate memory
    d->edgeCost = aligned_alloc(32, d->nodesCount * d->nodesCount * sizeof(double) + 4 * sizeof(double));

    // init data structure to pass to threads
    ThreadedInstance thInst = { .d = d, .freeElement = 0 };
    pthread_mutex_init (&thInst.mutex, NULL);

    // define thWorspace (number of matrix elements that a thread "locks" while computing distance on them)
    if (d->nodesCount > MAX_THWORSKPACE)   // very big matrix, lock fixed number
        thInst.thWorkspace = MAX_THWORSKPACE / d->params.threadsCount;
    else    // if there are just a few lines each thread "locks" from the begining <N°_OF_LINES>/<N°_OF_THREADS>
        thInst.thWorkspace = d->nodesCount / d->params.threadsCount;
    thInst.thWorkspace += 4 - (thInst.thWorkspace % 4); // make a multiple of 4 -> can use aligned instructions
    LOG (LOG_LVL_NOTICE, "Distance Matrix Computation: Thread Workspace value = %ld", thInst.thWorkspace);

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        pthread_create (threads[i], NULL, computeSqDistMatThread, &thInst);
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

    size_t thWorkspaceBegin = 0;
    size_t n = th->d->nodesCount;
    int finish = 0;

    register __m256d x1, y1, x2, y2, dist;
    register __m256i mask;

    while (finish == 0)
    {
        // CRITICAL REGION BEGIN
        // here thread gets it's workspace
        pthread_mutex_lock (&th->mutex);

        thWorkspaceBegin = th->freeElement;
        th->freeElement += th->thWorkspace;

        pthread_mutex_unlock (&th->mutex);
        // CRITICAL REGION END

        // now the thread can compute the distance matrix inside it's workspace
        
        

    }
    
    
}
    