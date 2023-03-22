#include "tsp.h"

typedef struct 
{
    Instance *d;
    //pthread_mutex_t mutex;
    pthread_mutex_t mutexNextEdgeMutex;
    pthread_mutex_t mutexBestConfigUpdate;
    pthread_mutex_t mutexFinishFlagMutex;
    size_t nextEdge;
    char finishedFlag;
    int *expandedBestSol; // copy of d.solution.bestSol with an extra element at the end that is bestSol[0] to close the tour
    float *bestCostOffset; // best optimization cost offset stored here. When negative a cost-reducing optimization has been found otherwise it's positive(or a very small negative to prevent numerical issues)
    size_t bestOffsetIDs[2]; // best optimization configuration: ids of the nodes that give the best costOffset when "switched" with 2opt
} ThreadedInstance;


static double _2optBestFix(Instance *d);

static void * _2optBestFixThread(void * arg);

static double _2optFirstFound(Instance *d);


static double _2optBestFix(Instance *d)
{
    ThreadedInstance th = { .d = d, .nextEdge = 0, .finishedFlag = 0 };

    // init expandedBestSol
    th.expandedBestSol = malloc((d->nodesCount + 1) * sizeof(int));
    memcpy(th.expandedBestSol, d->solution.bestSolution, d->nodesCount * sizeof(int));
    th.expandedBestSol[d->nodesCount] = d->solution.bestSolution[0];

    // init mutex
    pthread_mutex_init(&th.mutexFinishFlagMutex, NULL);
    pthread_mutex_init(&th.mutexNextEdgeMutex, NULL);

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
    size_t n = th->d->nodesCount;

    while ((pthread_mutex_lock(&th->mutexFinishFlagMutex) == 0) && (th->finishedFlag == 0)) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        pthread_mutex_unlock(&th->mutexFinishFlagMutex);
        // CRITICAL REGION for the finishedFlag: starts in the while condition and ends here

        // setup local values to avoid calling mutex too often
        float localBestCostOffset = INFINITY;
        size_t localBestOffsetID[2];

        // do while there are edges to be checked
        while ((pthread_mutex_lock(&th->mutexNextEdgeMutex)) && (th->nextEdge < n)) // check for one edge at a time every other edge(except already checked)
        {
            size_t edgeID = th->nextEdge; // goes from 0 (means edge from first node in solution to second) to n which is the number of nodes
            th->nextEdge++;
            pthread_mutex_unlock(&th->mutexNextEdgeMutex);
            // CRITICAL REGION for nextLine value: starts in the nested while condition end ends here

            float solEdgeWgt1 = th->d->edgeCost.mat[th->expandedBestSol[edgeID] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[edgeID + 1]];
            for (size_t i = 1 + edgeID; i < n; i++)
            {
                float solEdgeWgt2 = th->d->edgeCost.mat[th->expandedBestSol[i] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[i+1]];

                // check the combined weight other combination of edges
                float altEdge1 = th->d->edgeCost.mat[th->expandedBestSol[edgeID] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[i+1]];
                float altEdge2 = th->d->edgeCost.mat[th->expandedBestSol[i] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[edgeID+1]];
                if (solEdgeWgt1 + solEdgeWgt2 > altEdge1 + altEdge2)
                {
                    // update local best if necessary
                    if (localBestCostOffset < altEdge1 + altEdge2)
                    {
                        localBestCostOffset = (altEdge1 + altEdge2) - (solEdgeWgt1 + solEdgeWgt2);
                        localBestOffsetID[0] = edgeID;
                        localBestOffsetID[1] = i;
                    }
                }
            }
            
        }
        // update solution with best occurence found -> only one thread gets here, the others wait for mutex and then find the nextEdge = 0

        pthread_mutex_unlock(&th->mutexNextEdgeMutex);
    }
    



    pthread_mutex_unlock(&th->mutexFinishFlagMutex);
}