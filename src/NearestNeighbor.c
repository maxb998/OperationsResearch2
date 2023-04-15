#include "NearestNeighbor.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


typedef struct
{
    Solution *sol;
    pthread_mutex_t nodeLock;
    pthread_mutex_t saveLock;
    size_t startingNode;
}ThreadsData;

static void * threadNN(void *thInst);

// finds the closest unvisited node (pathCost is also updated in this method)
static inline size_t findSuccessorID(float *X, float *Y, size_t iterNode, Instance *inst, double *pathCost);

static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2);

Solution NearestNeighbor(Instance *inst)
{
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    Solution sol = newSolution(inst);

    // we initialize the seed if it has been passed as argument
    if(inst->params.randomSeed != -1) srand(inst->params.randomSeed);

    // we create and initialize the threaded instance
    ThreadsData th = {.sol = &sol, .startingNode = 0};

    pthread_mutex_init(&th.nodeLock, NULL);
    pthread_mutex_init(&th.saveLock, NULL);
    
    pthread_t threads[inst->params.nThreads];
    for(int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_create(&threads[i], NULL, threadNN, &th);
        LOG(LOG_LVL_DEBUG, "Nearest Neighbour : Thread %d CREATED", i);
    }

    for(int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join(threads[i], NULL);
        LOG(LOG_LVL_DEBUG, "Nearest Neighbour : Thread %d finished", i);
    }

    pthread_mutex_destroy(&th.nodeLock);
    pthread_mutex_destroy(&th.saveLock);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    sol.execTime = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0);

    return sol;
}

static void * threadNN(void *thInst)
{
    ThreadsData *th = (ThreadsData *)thInst;
    Solution *sol = th->sol;
    Instance *inst = sol->instance;

    // Allocate memory to contain the work-in-progress solution
    Solution tempSol = newSolution(inst);

    // We want the threads to repeat the computation for every node
    // For this we use a mutex on startingNode until it reaches nNodes
    while((pthread_mutex_lock(&th->nodeLock) == 0) && (th->startingNode < inst->nNodes))
    {
        // name of the starting node for this iteration of NearestNeighbour
        size_t iterNode = th->startingNode;

        th->startingNode++;

        // after incrementing the value of startingNode we can unlock the mutex
        pthread_mutex_unlock(&th->nodeLock);

        // copy all the nodes in the original order from instance into tempSol.X/Y (Take advantage of the way the memory is allocated for the best and current sol)
        for (size_t i = 0; i < (inst->nNodes + AVX_VEC_SIZE) * 2; i++)
            tempSol.X[i] = inst->X[i];

        // reset tempSol.indexPath to match the original
        for (int i = 0; i < (int)inst->nNodes + 1; i++)
            tempSol.indexPath[i] = i;

        // set first element of tempSol.X/Y to the element at index startingNode -> swap pos 0 with starting node
        swapElementsInSolution(&tempSol, 0, iterNode);

        // initialize the cost of the path to zero
        double pathCost = 0.;

        for(size_t i = 1; i < inst->nNodes - 1; i++)    // for n nodes we want to run this loop n-1 times, at the end we set as successor of the last node the starting node
        {
            size_t successorID = findSuccessorID(tempSol.X, tempSol.Y, i, inst, &pathCost);

            // Control on validity of successor: must be in [0,nNodes)
            if (successorID < 0 || successorID >= inst->nNodes)
                throwError(NULL, sol, "threadNN: error computing successor %d, value returned: %d", i, successorID);

            // set successor by swapping the element corresponding to succesorID with element i
            if (successorID != i)
            {
                {   // swap coordinates
                    register float temp;
                    swapElems(tempSol.X[i], tempSol.X[successorID], temp);
                    swapElems(tempSol.Y[i], tempSol.Y[successorID], temp);
                }

                { // swap index
                    register int temp;
                    swapElems(tempSol.indexPath[i], tempSol.indexPath[successorID], temp);
                }
            }
        }
        // at the end we set the starting node as successor of the last one to close the circuit
        tempSol.X[inst->nNodes] = tempSol.X[0];
        tempSol.Y[inst->nNodes] = tempSol.Y[0];
        tempSol.indexPath[inst->nNodes] = tempSol.indexPath[0];

        // add last and previous to last edge weights
        pathCost += squaredEdgeCost(tempSol.X[inst->nNodes-2], tempSol.Y[inst->nNodes-2], tempSol.X[inst->nNodes-1], tempSol.Y[inst->nNodes-1], inst->params.edgeWeightType);
        pathCost += squaredEdgeCost(tempSol.X[inst->nNodes-1], tempSol.Y[inst->nNodes-1], tempSol.X[inst->nNodes], tempSol.Y[inst->nNodes], inst->params.edgeWeightType);

        // to check if we have to update the best solution we use another mutex
        if((pthread_mutex_lock(&th->saveLock) == 0) && (sol->bestCost > pathCost))
        {
            sol->bestCost = pathCost;

            { // swap pointers with a macro(swapPtr) to declutter code (macro is defined at the top)
                register float *temp;
                swapElems(sol->X, tempSol.X, temp);
                swapElems(sol->Y, tempSol.Y, temp);
            }
            {
                register int *temp;
                swapElems(sol->indexPath, tempSol.indexPath, temp);
            }

            LOG(LOG_LVL_EVERYTHING, "Found better solution starting from node %ld, cost: %f", iterNode, pathCost);
        }
        pthread_mutex_unlock(&th->saveLock);
    }

    pthread_mutex_unlock(&th->nodeLock);

    destroySolution(&tempSol);

    pthread_exit(NULL);
}

static inline size_t findSuccessorID(float *X, float *Y, size_t iterNum, Instance *inst, double *pathCost)
{
    float minVecStore[AVX_VEC_SIZE];
    int IDVecStore[AVX_VEC_SIZE];

    // to keep track of the first best distance
    __m256 min1Vec = _mm256_set1_ps(INFINITY);
    __m256i min1IDsVec = _mm256_set1_epi32(-1);

    // to keep track of the second best distance
    //    (if the second best distance gets loaded in the same position in the avx vector as the best distance,
    //     otherwise the second best distance is still inside the vectors defined above (P(2ndBest in above vec) = 1 - 1 / AVX_VEC_SIZE))
    //__m256 min2Vec = _mm256_set1_ps(INFINITY);
    //__m256i min2IDsVec = _mm256_set1_epi32(-1);

    //__m256 infinite = _mm256_set1_ps(INFINITY);

    __m256i ids = _mm256_setr_epi32(0,1,2,3,4,5,6,7), increment = _mm256_set1_epi32(AVX_VEC_SIZE);
    __m256 x1 = _mm256_set1_ps(X[iterNum-1]), y1 = _mm256_set1_ps(Y[iterNum-1]);
    
    // set ids to right starting position
    ids = _mm256_add_epi32(ids, _mm256_set1_epi32((int)iterNum));

    // check if we are working with rounded weights
    for (size_t i = iterNum; i < inst->nNodes; i += AVX_VEC_SIZE)
    {
        // get distance first
        __m256 x2 = _mm256_loadu_ps(&X[i]), y2 = _mm256_loadu_ps(&Y[i]);
        __m256 dist = squaredEdgeCost_VEC(x1, y1, x2, y2, inst->params.edgeWeightType);

        // get first minimum
        __m256i mask = _mm256_castps_si256(_mm256_cmp_ps(dist, min1Vec, _CMP_LT_OQ)); // set for "non-signaling" -> if one element is NaN than set as false
        min1IDsVec = _mm256_blendv_epi8(min1IDsVec, ids, mask);  // get 1stMin id
        min1Vec = _mm256_min_ps(min1Vec, dist); // get 1st minumum value

        // replace elements that has been inserted into min1Vec with INFINITE so that they don't get inserted in min2Vec as well
        //dist = _mm256_blendv_ps(dist, infinite, _mm256_castsi256_ps(mask));

        // get second minimum
        //mask = _mm256_castps_si256(_mm256_cmp_ps(dist, min2Vec, _CMP_LT_OQ));
        //min2IDsVec = _mm256_blendv_epi8(min2IDsVec, ids, mask);
        //min2Vec = _mm256_min_ps(min2Vec, dist);

        // increment ids
        ids = _mm256_add_epi32(ids, increment);
    }

    // store vectors
    _mm256_storeu_ps(minVecStore, min1Vec);
    //_mm256_store_ps(&minVecStore[AVX_VEC_SIZE], min2Vec);
    _mm256_storeu_si256((__m256i*)IDVecStore, min1IDsVec);
    //_mm256_store_si256((__m256i*)&IDVecStore[AVX_VEC_SIZE], min2IDsVec);

    // find id of minimum and second minimum in 8 elements
    size_t minIDID = -1, secondMinIDID = -1; // are IDs of that must be placed int IDVecStore[] to get the real minID so they are IDs of IDs(? not sure if explanation is good but it works)
    float min = INFINITY, secondMin = INFINITY;
    for (size_t i = 0; i < AVX_VEC_SIZE /* * 2 */; i++)
    {
        if (min > minVecStore[i])
        {
            secondMin = min;
            secondMinIDID = minIDID;
            min = minVecStore[i];
            minIDID = i;
        }
        else if (secondMin > minVecStore[i])
        {
            secondMin = minVecStore[i];
            secondMinIDID = i;
        }
    }

    // We choose what node of the two best we return if GRASP has been required
    if (inst->params.randomSeed != -1 && rand() > NN_GRASP_COEFF && IDVecStore[secondMinIDID] != -1)
    {
        *pathCost += minVecStore[secondMinIDID];
        return (size_t)IDVecStore[secondMinIDID];
    }
    else
    {
        *pathCost += minVecStore[minIDID];
        return (size_t)IDVecStore[minIDID];
    }
}


static inline void swapElementsInSolution(Solution *sol, size_t pos1, size_t pos2)
{
    register float tempf;
    swapElems(sol->X[pos1], sol->X[pos2], tempf);
    swapElems(sol->Y[pos1], sol->Y[pos2], tempf);
    register int tempi;
    swapElems(sol->indexPath[pos1], sol->indexPath[pos2], tempi);
}