//#include "tsp.h"
#include "edgeCostFunctions.h"


typedef struct
{
    Solution *sol;
    pthread_mutex_t nodeLock;
    pthread_mutex_t saveLock;
    size_t startingNode;
}ThreadedInstance;

static void * threadNN(void *thInst);

// finds the closest unvisited node (pathCost is also updated in this method)
static inline size_t findSuccessorID(float *X, float *Y, size_t iterNode, Instance *inst, float *pathCost, float *minVecStore, int *IDVecStore);

static inline void swapf(float *arr, size_t pos1, size_t pos2);

Solution NearestNeighbour(Instance *inst)
{
    clock_t start, end;
    start = clock();

    Solution sol = newSolution(inst);

    // we initialize the seed if it has been passed as argument
    if(inst->params.randomSeed != -1) srand(inst->params.randomSeed);

    // we create and initialize the threaded instance
    ThreadedInstance th = {.sol = &sol, .startingNode = 0};

    pthread_mutex_init(&th.nodeLock, NULL);
    pthread_mutex_init(&th.saveLock, NULL);
    
    pthread_t threads[inst->params.nThreads];
    for(int i = 0; i < inst->params.nThreads; i++)
    {
        pthread_create(&threads[i], NULL, threadNN, &th);
        LOG(LOG_LVL_LOG, "Nearest Neighbour : Thread %d CREATED", i);
    }

    double maxThreadTime = 0.0;
    for(int i = 0; i < inst->params.nThreads; i++)
    {
        double *returnPtr;
        pthread_join(threads[i], (void **)&returnPtr);

        double threadTime = *returnPtr;
        free(returnPtr);

        if (maxThreadTime < threadTime)
            maxThreadTime = threadTime;

        LOG(LOG_LVL_LOG, "Nearest Neighbour : Thread %d finished in %e time", i, threadTime);
    }

    pthread_mutex_destroy(&th.nodeLock);
    pthread_mutex_destroy(&th.saveLock);

    end = clock();
    sol.execTime = maxThreadTime + ((double)(end - start))/CLOCKS_PER_SEC;

    return sol;
}

static void * threadNN(void *thInst)
{
    clock_t start, end;
    start = clock();

    ThreadedInstance *th = (ThreadedInstance *)thInst;
    Solution *sol = th->sol;
    Instance *inst = sol->instance;

    // Allocate memory to contain the work-in-progress solution
    float *currentSolX = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float));
    float *currentSolY = &currentSolX[inst->nNodes + AVX_VEC_SIZE];

    // Allocate memory to store Vector register element (aligned for not apparent reason than it feels better)
    float *minVecStore = aligned_alloc(32, 2 * AVX_VEC_SIZE * sizeof(float));
    int *IDVecStore = aligned_alloc(32, 2 * AVX_VEC_SIZE * sizeof(int));

    // We want the threads to repeat the computation for every node
    // For this we use a mutex on startingNode until it reaches nNodes
    while((pthread_mutex_lock(&th->nodeLock) == 0) && (th->startingNode < inst->nNodes))
    {
        // name of the starting node for this iteration of NearestNeighbour
        size_t iterNode = th->startingNode;

        th->startingNode++;

        // after incrementing the value of startingNode we can unlock the mutex
        pthread_mutex_unlock(&th->nodeLock);


        // copy all the nodes in the original order from instance into currentSolX/Y (Take advantage of the way the memory is allocated for the best and current sol)
        memcpy(currentSolX, inst->X, (inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float)); // this also copies Y

        // set first element of currentSolX/Y to the element at index startingNode -> swap pos 0 with starting node
        swapf(currentSolX, 0, iterNode);
        swapf(currentSolY, 0, iterNode);

        // reset last element
        currentSolX[inst->nNodes] = INFINITY;
        currentSolY[inst->nNodes] = INFINITY;

        // initialize the cost of the path to zero
        float pathCost = 0.;

        for(size_t i = 1; i < inst->nNodes - 1; i++)    // for n nodes we want to run this loop n-1 times, at the end we set as successor of the last node the starting node
        {
            size_t successorID = findSuccessorID(currentSolX, currentSolY, i, inst, &pathCost, minVecStore, IDVecStore);

            // Control on validity of successor: must be in [0,nNodes)
            if (successorID < 0 || successorID >= inst->nNodes)
                throwError(NULL, sol, "threadNN: error computing successor %d, value returned: %d", i, successorID);

            // set successor by swapping the element corresponding to succesorID with element i
            if (successorID != i)
            {
                swapf(currentSolX, i, successorID);
                swapf(currentSolY, i, successorID);
            }
        }
        // at the end we set the starting node as successor of the last one to close the circuit
        currentSolX[inst->nNodes] = currentSolX[0];
        currentSolY[inst->nNodes] = currentSolY[0];

        if (th->startingNode == 41)
            printf("HERE\n");

        // add last and previours to last edge weights
        pathCost += squaredEdgeCost(currentSolX[inst->nNodes-2], currentSolY[inst->nNodes-2], currentSolX[inst->nNodes-1], currentSolY[inst->nNodes-1], inst->params.edgeWeightType);
        pathCost += squaredEdgeCost(currentSolX[inst->nNodes-1], currentSolY[inst->nNodes-1], currentSolX[inst->nNodes], currentSolY[inst->nNodes], inst->params.edgeWeightType);

        // to check if we have to update the best solution we use another mutex
        if((pthread_mutex_lock(&th->saveLock) == 0) && (sol->bestCost > pathCost))
        {
            float rnd = squaredEdgeCost(currentSolX[inst->nNodes-1], currentSolY[inst->nNodes-1], currentSolX[inst->nNodes], currentSolY[inst->nNodes], inst->params.edgeWeightType);
            printf("%f\n", rnd);

            sol->bestCost = pathCost;

            float * swapPtr = sol->X;
            sol->X = currentSolX;
            currentSolX = swapPtr;

            swapPtr = sol->Y;
            sol->Y = currentSolY;
            currentSolY = swapPtr;

            LOG(LOG_LVL_LOG, "Found better solution starting from node %ld, cost: %f", iterNode, pathCost);
        }
        pthread_mutex_unlock(&th->saveLock);
    }

    pthread_mutex_unlock(&th->nodeLock);


    free(minVecStore);
    free(IDVecStore);

    double *threadTime = malloc(sizeof(double));
    end = clock();
    *threadTime = ((double)(end - start))/CLOCKS_PER_SEC;

    pthread_exit((void*)threadTime);
}

static inline size_t findSuccessorID(float *X, float *Y, size_t iterNum, Instance *inst, float *pathCost, float *minVecStore, int *IDVecStore)
{
    // to keep track of the first best distance
    __m256 min1Vec = _mm256_set1_ps(INFINITY);
    __m256i min1IDsVec = _mm256_set1_epi32(-1);

    // to keep track of the second best distance
    //    (if the second best distance gets loaded in the same position in the avx vector as the best distance,
    //     otherwise the second best distance is still inside the vectors defined above (P(2ndBest in above vec) = 1 - 1 / AVX_VEC_SIZE))
    //__m256 min2Vec = _mm256_set1_ps(INFINITY);
    //__m256i min2IDsVec = _mm256_set1_epi32(-1);

    //__m256 infinite = _mm256_set1_ps(INFINITY);

    __m256i idsVec = _mm256_setr_epi32(0,1,2,3,4,5,6,7), increment = _mm256_set1_epi32(AVX_VEC_SIZE);
    __m256 x1Vec = _mm256_set1_ps(X[iterNum-1]), y1Vec = _mm256_set1_ps(Y[iterNum-1]);
    
    // set idsVec to right starting position
    idsVec = _mm256_add_epi32(idsVec, _mm256_set1_epi32((int)iterNum));

    // check if we are working with rounded weights
    for (size_t i = iterNum; i < inst->nNodes; i += AVX_VEC_SIZE)
    {
        // get distance first
        __m256 x2Vec = _mm256_loadu_ps(&X[i]), y2Vec = _mm256_loadu_ps(&Y[i]);
        __m256 dist = squaredEdgeCost_VEC(x1Vec, y1Vec, x2Vec, y2Vec, inst->params.edgeWeightType);

        // get first minimum
        __m256i mask = _mm256_castps_si256(_mm256_cmp_ps(dist, min1Vec, _CMP_LT_OQ)); // set for "non-signaling" -> if one element is NaN than set as false
        min1IDsVec = _mm256_blendv_epi8(min1IDsVec, idsVec, mask);  // get 1stMin id
        min1Vec = _mm256_min_ps(min1Vec, dist); // get 1st minumum value

        // replace elements that has been inserted into min1Vec with INFINITE so that they don't get inserted in min2Vec as well
        //dist = _mm256_blendv_ps(dist, infinite, _mm256_castsi256_ps(mask));

        // get second minimum
        //mask = _mm256_castps_si256(_mm256_cmp_ps(dist, min2Vec, _CMP_LT_OQ));
        //min2IDsVec = _mm256_blendv_epi8(min2IDsVec, idsVec, mask);
        //min2Vec = _mm256_min_ps(min2Vec, dist);

        // increment ids
        idsVec = _mm256_add_epi32(idsVec, increment);
    }

    // store vectors
    _mm256_store_ps(minVecStore, min1Vec);
    //_mm256_store_ps(&minVecStore[AVX_VEC_SIZE], min2Vec);
    _mm256_store_si256((__m256i*)IDVecStore, min1IDsVec);
    //_mm256_store_si256((__m256i*)&IDVecStore[AVX_VEC_SIZE], min2IDsVec);

    // find id of minimum and second minimum in 16 elements
    size_t minIDID = 0, secondMinIDID = 0; // are IDs of that must be placed int IDVecStore[] to get the real minID so they are IDs of IDs(? not sure if explanation is good but it works)
    float min = INFINITY;// secondMin = INFINITY;
    for (size_t i = 0; i < AVX_VEC_SIZE /* * 2 */; i++)
    {
        if (min > minVecStore[i])
        {
            //secondMin = min;
            //secondMinIDID = minIDID;
            min = minVecStore[i];
            minIDID = i;
        }
        /*else if (secondMin > minVecStore[i]) // also check that it's not the same value as minIDID
        {
            secondMin = minVecStore[i];
            secondMinIDID = i;
        }*/
    }

    // We choose what node of the two best we return if GRASP has been required
    if (inst->params.randomSeed != -1 && rand() > GRASP_COEFF && IDVecStore[secondMinIDID] != -1)
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

static inline void swapf(float *arr, size_t pos1, size_t pos2)
{
    float swapVal = arr[pos1];
    arr[pos1] = arr[pos2];
    arr[pos2] = swapVal;
}