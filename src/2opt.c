#include "tsp.h"

typedef struct 
{
    Instance *d;
    //pthread_mutex_t mutex;
    pthread_mutex_t mutexNextLineMutex;
    pthread_mutex_t mutexFinishFlagMutex;
    size_t nextPtID;
    char finishedFlag;

} ThreadedInstance;


static double _2optBestFix(Instance *d);

static void * _2optBestFixThread(void * arg);

static double _2optFirstFound(Instance *d);


static double _2optBestFix(Instance *d)
{
    ThreadedInstance th = { .d = d, .nextPtID = 0, .finishedFlag = 0 };

    // init mutex
    pthread_mutex_init(&th.mutexFinishFlagMutex, NULL);
    pthread_mutex_init(&th.mutexNextLineMutex, NULL);

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < d->params.threadsCount; i++)
        pthread_create(&threads[i], NULL, _2optBestFixThread, &th);

    // wait for threads and get execution time
    double threadsTime[MAX_THREADS];
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        double * returnedTime;
        pthread_join(threads[i], &returnedTime);

        if (!returnedTime) throwError(d, "Return value of thread %ld int 2optBestFix does not exist", i);

        threadsTime[i] = *returnedTime;
        free(returnedTime);
    }
    
    
}

static void * _2optBestFixThread(void * arg)
{
    ThreadedInstance *th = (ThreadedInstance*)arg;

    while ((pthread_mutex_lock(&th->mutexFinishFlagMutex) == 0) && (th->finishedFlag == 0))
    {
        pthread_mutex_unlock(&th->mutexFinishFlagMutex);
        
    }
    pthread_mutex_unlock(&th->mutexFinishFlagMutex);
}