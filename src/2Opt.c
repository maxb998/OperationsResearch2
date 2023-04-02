#include "2Opt.h"
#include "TspUtilities.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time



typedef struct
{
    Solution *sol;

    // id of the next edge to be checked against all solution edges
    size_t nextEdge;

    // flag that is set to 1 when the bestSolution has stopped improving
    char finishedFlag;

    // best optimization cost offset stored here.
    // When negative a cost-reducing optimization has been found otherwise it's positive(or a very small negative to prevent numerical issues)
    float bestOffset;

    // best optimization configuration: ids of the nodes that give the best costOffset when "switched" with 2opt
    size_t bestOffsetIDs[2];

    // Each thread set threadFinishFlag[threadID] = 1 when it's waiting for other thread to finish and update the bestSolution
    int threadFinishFlag[MAX_THREADS];

    // flag to use in combination with conditionBestSolUpdate to wait for the thread currently updating the best solution
    //int bestSolUpdateFinished;

    // mutex
    pthread_mutex_t mutex;

    // condition
    pthread_cond_t conditionBestSolUpdated;

} ThreadsData;

// struct that allows to give a meaningful and unique id to each thread
typedef struct 
{
    ThreadsData *th;
    int threadID;
} ThreadsDataWithID;


static inline void bestSolutionUpdate(ThreadsData *th);

static void *_2optBestFixCostMatrixThread(void *arg);



double _2optBestFix(Solution *sol)
{
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    ThreadsData th = { .sol = sol, .nextEdge = 0, .finishedFlag = 0, .bestOffset = 0 };
    Instance *inst = sol->instance;

    // copy first element to last position (memory allocation is AVX_VEC_SIZE elemens loose)
    sol->indexPath[inst->nNodes] = sol->indexPath[0];

    // init mutexes & condition
    pthread_mutex_init(&th.mutex, NULL);
    pthread_cond_init(&th.conditionBestSolUpdated, NULL);

    // struct necessary to give each thread it's identifier correctly
    ThreadsDataWithID thIDs[MAX_THREADS];
    for (size_t i = 0; i < inst->params.nThreads; i++)
    {
        thIDs[i].th = &th;
        thIDs[i].threadID = i;
    }

    // start threads
    pthread_t workers[MAX_THREADS];
    for (size_t i = 0; i < inst->params.nThreads; i++)
    {
        pthread_create(&workers[i], NULL, _2optBestFixCostMatrixThread, &thIDs[i]);
        LOG(LOG_LVL_DEBUG, "2-Opt: Thread %ld created", i);
    }

    // wait for threads
    for (size_t i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join(workers[i], NULL);
        LOG(LOG_LVL_DEBUG, "2-Opt: Thread %ld finished", i);
    }

    // destroy mutexes and conditions
    pthread_mutex_destroy(&th.mutex);
    pthread_cond_destroy(&th.conditionBestSolUpdated);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    double elapsed = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);
    return elapsed;
}

static void *_2optBestFixCostMatrixThread(void *arg)
{
    ThreadsDataWithID *thID = (ThreadsDataWithID*)arg;
    ThreadsData *th = thID->th;
    int threadID = thID->threadID;
    Solution *sol = th->sol;
    Instance *inst = sol->instance;

    // setup "shortcuts" variables to declutter the code
    size_t n = inst->nNodes;

    while (th->finishedFlag == 0) // runs 2opt until no more moves are made in one iteration of 2opt
    {
        // setup local values to avoid calling mutex too often
        float localBestOffset = 0;
        size_t localBestOffsetID[2];

        // CRITICAL SECTION START
        // do while there are edges to be checked
        while ((pthread_mutex_lock(&th->mutex) == 0) && (th->nextEdge < n - 1)) // check for one edge at a time every other edge(except already checked)
        {
            size_t edgeID = th->nextEdge; // goes from 0 (means edge from first node in solution to second) to n which is the number of nodes
            th->nextEdge++;
            pthread_mutex_unlock(&th->mutex);
            // // CRITICAL SECTION END

            float solEdgeWgt1 = inst->edgeCostMat[sol->indexPath[edgeID] * n + sol->indexPath[edgeID + 1]];
            for (size_t i = 2 + edgeID; (i < n - 1) || ((i < n) && (edgeID > 0)); i++)
            {
                float solEdgeWgt2 = inst->edgeCostMat[sol->indexPath[i] * n + sol->indexPath[i + 1]];

                // check the combined weight other combination of edges
                float altEdge1 = inst->edgeCostMat[sol->indexPath[edgeID] * n + sol->indexPath[i]];
                float altEdge2 = inst->edgeCostMat[sol->indexPath[i + 1] * n + sol->indexPath[edgeID + 1]];
                if (solEdgeWgt1 + solEdgeWgt2 > altEdge1 + altEdge2)
                {
                    // update local best if current one is better
                    float currentOffset = (altEdge1 + altEdge2) - (solEdgeWgt1 + solEdgeWgt2);
                    if (localBestOffset > currentOffset)
                    {
                        localBestOffset = currentOffset;
                        localBestOffsetID[0] = edgeID;
                        localBestOffsetID[1] = i;
                    }
                }
            }
        }
        // When this point is reached the mutex is still locked: STILL IN CRITICAL SECTION

        // When all possibilities for one edge have been checked, update thread shared bestoffset if localbest is better
        if (th->bestOffset > localBestOffset)
        {
            th->bestOffset = localBestOffset;
            th->bestOffsetIDs[0] = localBestOffsetID[0];
            th->bestOffsetIDs[1] = localBestOffsetID[1];
        }

        th->nextEdge = 0;
        th->threadFinishFlag[threadID] = 1;

        int lastThreadFlag = 1; // if it remains one than this is the last thread going through this part at the end of the 2opt move search
        for (size_t i = 0; i < inst->params.nThreads; i++)
            if (th->threadFinishFlag[i] == 0)
                lastThreadFlag = 0;
        

        // signal updater thread to update solution only when all worker threads have finished
        if (lastThreadFlag == 1)
        {
            // update the best solution
            bestSolutionUpdate(th);

            // reset all the threadFinishFlag to 0
            for (size_t i = 0; i < inst->params.nThreads; i++)
                th->threadFinishFlag[i] = 0;

            // signal other threads that bestSolution update has been done
            pthread_cond_broadcast(&th->conditionBestSolUpdated);
        }
        else
            // if not the last thread to finish wait
            while (th->threadFinishFlag[threadID] == 1)
                pthread_cond_wait(&th->conditionBestSolUpdated, &th->mutex);

        pthread_mutex_unlock(&th->mutex);
        // CRITICAL SECTION END
    }

    pthread_exit(NULL);
}


static inline void bestSolutionUpdate(ThreadsData *th)
{
    Solution *sol = th->sol;
    if (th->bestOffset < -EPSILON) // avoid to do practically meaningless optimization and float precision errors(TO REVISE AND DO A RELATIVE COMPARISON)
    {
        /*
         * UPDATE BEST SOLUTION ########################################################################################
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
            int temp = sol->indexPath[smallID];
            sol->indexPath[smallID] = sol->indexPath[bigID];
            sol->indexPath[bigID] = temp;
            smallID++;
            bigID--;
        }

        th->bestOffset = 0;
    }
    else // NO COST IMPROVING 2opt MOVE HAS BEEN FOUND
        th->finishedFlag = 1;
}
