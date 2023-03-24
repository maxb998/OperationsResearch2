#include "tsp.h"

typedef struct 
{
    Instance *d;

    // id of the next edge to be checked against all solution edges
    size_t nextEdge;

    // flag that is set to 1 when the bestSolution has stopped improving
    char finishedFlag;

    // copy of d.solution.bestSol with an extra element at the end that is bestSol[0] to close the tour
    int *expandedBestSol;

    // best optimization cost offset stored here.
    // When negative a cost-reducing optimization has been found otherwise it's positive(or a very small negative to prevent numerical issues)
    float bestCostOffset; 

    // best optimization configuration: ids of the nodes that give the best costOffset when "switched" with 2opt
    size_t bestOffsetIDs[2];

    // when each thread reach the end of one iteration to find the best edges to which apply 2opt gives this a + 1.
    // Needed when updating the bestSolution because we want to be 100% sure that all threads have finished the subroutine
    // Also needed to avoid condition release too much
    int threadFinish;

    // mutexes
    pthread_mutex_t mutexNextEdgeMutex;
    pthread_mutex_t mutexBestConfigUpdate;
    pthread_mutex_t mutexFinishFlagMutex;

    // condition
    pthread_cond_t conditionBestSolUpdate;

} ThreadedInstance;

static void * _2optBestFixThread(void * arg);


double _2optBestFix(Instance *d)
{
    ThreadedInstance th = { .d = d, .nextEdge = 0, .finishedFlag = 0, .threadFinish = 0 };

    // init expandedBestSol
    th.expandedBestSol = malloc((d->nodesCount + 1) * sizeof(int));
    memcpy(th.expandedBestSol, d->solution.bestSolution, d->nodesCount * sizeof(int));
    th.expandedBestSol[d->nodesCount] = d->solution.bestSolution[0];

    // init mutexes & condition
    pthread_mutex_init(&th.mutexFinishFlagMutex, NULL);
    pthread_mutex_init(&th.mutexBestConfigUpdate, NULL);
    pthread_mutex_init(&th.mutexNextEdgeMutex, NULL);
    pthread_cond_init(&th.conditionBestSolUpdate, NULL);

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < d->params.threadsCount; i++)
        pthread_create(&threads[i], NULL, _2optBestFixThread, &th);

    // wait for threads and get execution time
    double maxThreadsTime = 0.0;
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        double * returnedTime;
        pthread_join(threads[i], (void**)&returnedTime);

        if (!returnedTime) throwError(d, "Return value of thread %ld int 2optBestFix does not exist", i);

        if (maxThreadsTime < *returnedTime)
            maxThreadsTime = *returnedTime;
        free(returnedTime);
    }

    // destroy mutexes and conditions
    pthread_mutex_destroy(&th.mutexFinishFlagMutex);
    pthread_mutex_destroy(&th.mutexBestConfigUpdate);
    pthread_mutex_destroy(&th.mutexNextEdgeMutex);
    pthread_cond_destroy(&th.conditionBestSolUpdate);

    return 0.0;
}


static void * _2optBestFixThread(void * arg)
{
    clock_t start, end;
    start = clock();

    // CRITICAL SECTION #1 -> finishFlag read
    // CRITICAL SECTION #2 -> 
    // CRITICAL SECTION #3 -> 

    ThreadedInstance *th = (ThreadedInstance*)arg;
    size_t n = th->d->nodesCount;

    // CRITICAL SECTION #1 -> BEGIN (inside the while condition)
    while ((pthread_mutex_lock(&th->mutexFinishFlagMutex) == 0) && (th->finishedFlag == 0)) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        pthread_mutex_unlock(&th->mutexFinishFlagMutex);
        // CRITICAL SECTION #1 -> END # 1.1

        // CRITICAL SECTION #2 -> START
        // do while there are edges to be checked
        while ((pthread_mutex_lock(&th->mutexNextEdgeMutex)) && (th->nextEdge < n)) // check for one edge at a time every other edge(except already checked)
        {
            size_t edgeID = th->nextEdge; // goes from 0 (means edge from first node in solution to second) to n which is the number of nodes
            th->nextEdge++;
            pthread_mutex_unlock(&th->mutexNextEdgeMutex);
            // // CRITICAL SECTION #2 -> END # 2.1

            // setup local values to avoid calling mutex too often
            float localBestCostOffset = INFINITY;
            size_t localBestOffsetID[2];

            float solEdgeWgt1 = th->d->edgeCost.mat[th->expandedBestSol[edgeID] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[edgeID + 1]];
            for (size_t i = 1 + edgeID; i < n; i++)
            {
                float solEdgeWgt2 = th->d->edgeCost.mat[th->expandedBestSol[i] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[i+1]];

                // check the combined weight other combination of edges
                float altEdge1 = th->d->edgeCost.mat[th->expandedBestSol[edgeID] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[i+1]];
                float altEdge2 = th->d->edgeCost.mat[th->expandedBestSol[i] * th->d->edgeCost.rowSizeMem + th->expandedBestSol[edgeID+1]];
                if (solEdgeWgt1 + solEdgeWgt2 > altEdge1 + altEdge2)
                {
                    // update local best if current one is better
                    if (localBestCostOffset < altEdge1 + altEdge2)
                    {
                        localBestCostOffset = (altEdge1 + altEdge2) - (solEdgeWgt1 + solEdgeWgt2);
                        localBestOffsetID[0] = edgeID;
                        localBestOffsetID[1] = i;
                    }
                }
            }
            
            // CRITICAL SECTION #3 -> BEGIN #############
            // When all possibilities for one edge have been checked, update thread shared bestoffset if localbest is better
            pthread_mutex_lock(&th->mutexBestConfigUpdate);

            if (localBestCostOffset < th->bestCostOffset)
            {
                th->bestCostOffset = localBestCostOffset;
                th->bestOffsetIDs[0] = localBestOffsetID[0];
                th->bestOffsetIDs[1] = localBestOffsetID[1];
            }

            pthread_mutex_unlock(&th->mutexBestConfigUpdate);
            // CRITICAL SECTION #3 -> END ##################
        }
        // Update solution with best occurence found
        th->threadFinish++;
        if (th->threadFinish == th->d->params.threadsCount) // allow best solution update only if we are certain to have the best possible
        {
            th->threadFinish = 1;

            if (th->bestCostOffset > -EPSILON) // avoid to do practically meaningless optimization and float precision errors(TO REVISE AND DO A RELATIVE COMPARISON)
            {
                /*UPDATE BEST SOLUTION ########################################################################################
                 *      bestSolIDs = { 0 1 2 3 4 5 6 7 8 9 }
                 * if bestSolution = { 2 5 8 7 6 9 4 1 0 3 } and bestOffsetIDs = { 3 8 }
                 *          -> means edges to swapped are (7,6) and (0,3) with (7,0) and (0,6)
                 * at the end of the swap, the new bestSolution will be { 2 5 8 7 0 1 4 9 6 3 }     ^         ^
                 *     old bestSolution = { 2 5 8 7 6 9 4 1 0 3 }   with original indexes = { 0 1 2 3 4 5 6 7 8 9 }
                 *     new bestSolution = { 2 5 8 7 0 1 4 9 6 3 }   with original indexes = { 0 1 2 3 8 7 6 5 4 9 }
                 * 
                 * Which means that we must invert the elements of bestSolution from index 3(not inlcuded) to index 8(included)
                 */

                size_t smallID = th->bestOffsetIDs[0] + 1, bigID = th->bestOffsetIDs[1];

                while (smallID < bigID)
                {
                    int temp = th->d->solution.bestSolution[smallID];
                    th->d->solution.bestSolution[smallID] = th->d->solution.bestSolution[bigID];
                    th->d->solution.bestSolution[bigID] = temp;
                    smallID++;
                    bigID++;
                }
            }
            else // NO COST IMPROVING 2opt MOVE FOUND
            {
                th->finishedFlag = 1;
            }

            pthread_cond_broadcast(&th->conditionBestSolUpdate);
        }
        else
        {
            while (th->threadFinish != 1)
                pthread_cond_wait(&th->conditionBestSolUpdate, &th->mutexNextEdgeMutex);
            // once here the solution has been updated by another thread

        }

        pthread_mutex_unlock(&th->mutexNextEdgeMutex);
        // CRITICAL SECTION #2 -> END # 2.2
    }

    pthread_mutex_unlock(&th->mutexFinishFlagMutex);
    // CRITICAL SECTION #1 -> END # 1.2

    end = clock();

    double *cpuTimeUsed = malloc(sizeof(double));
    *cpuTimeUsed = ((double)(end-start)) / CLOCKS_PER_SEC;

    pthread_exit(cpuTimeUsed);
}