#include "tsp.h"

typedef struct
{
    Instance *inst;
    pthread_mutex_t mutex;
    size_t nextRow; // only value that needs critical region, used to select the row on which the thread will work on
} ThreadedInstance;

// method to compute edges weights with multiple threads using AVX2 instructions
static void * computeDistMatThread(void * arg);

void printDistanceMatrix(Instance *inst, int showEndRowPlaceholder)
{
    size_t rowSizeToPrint;
    if (showEndRowPlaceholder == 1)
        rowSizeToPrint = inst->edgeCost.rowElems;
    else
        rowSizeToPrint = inst->nNodes;

    LOG(LOG_LVL_LOG, "Printing the distance matrix");

    for (size_t row = 0; row < inst->nNodes; row++)
    {
        printf(" ");
        for (size_t col = 0; col < rowSizeToPrint; col++)
            printf("%.3e ", inst->edgeCost.data[row * inst->edgeCost.rowElems + col]);
        printf("\n");
    }
}



double computeDistanceMatrix(Instance *inst)
{
    double totCpuTime = 0.;
    clock_t start, end;

    start = clock();

    // first check that distance type is supported
    if (inst->params.edgeWeightType > ATT)
        throwError(inst, NULL, "Distance Matrix Computation: Edge weight type unsupported");

    // allocate memory
    size_t n = inst->nNodes;

    if (n % AVX_VEC_SIZE > 0)
        inst->edgeCost.rowElems = n + AVX_VEC_SIZE - n % AVX_VEC_SIZE;   // way easier to work with avx instructions if each row start pointer is alligned
    else
        inst->edgeCost.rowElems = n;

    inst->edgeCost.data = malloc(inst->edgeCost.rowElems * n * sizeof(float));
        
    // init data structure to pass to threads
    ThreadedInstance thInst = { .inst = inst, .nextRow = 0 };
    pthread_mutex_init (&thInst.mutex, NULL);

    // each thread compute one row of the distance matrix before calling mutex to get a new workspace

    // start threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < inst->params.nThreads; i++)
    {
        pthread_create (&threads[i], NULL, computeDistMatThread, &thInst);
        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld CREATED", i);
    }
    
    end = clock();
    totCpuTime += ((double)(end-start)) / CLOCKS_PER_SEC;

    // wait threads to finish and get time of each thread
    double threadsTime[MAX_THREADS], *returnPtr;
    for (size_t i = 0; i < inst->params.nThreads; i++)
    {
        pthread_join (threads[i], (void**)&returnPtr);
        threadsTime[i] = *returnPtr;
        free(returnPtr);

        LOG (LOG_LVL_DEBUG, "Distance Matrix Computation: Thread %ld finished in %.3e seconds", i, threadsTime[i]);
    }

    double threadMax = 0.;
    for (size_t i = 0; i < inst->params.nThreads; i++)
        if (threadMax < threadsTime[i])
            threadMax = threadsTime[i];

    return totCpuTime + threadMax;
}

// EUCLIDEAN DISTANCE
static void * computeDistMatThread(void* d)
{
    // measure first tick
    clock_t start, end;
    start = clock();

    // correctly identify data
    ThreadedInstance *th = (ThreadedInstance*)d;

    size_t n = th->inst->nNodes;

    while ( (pthread_mutex_lock(&th->mutex) == 0) && (th->nextRow < th->inst->nNodes) )    // lock mutex before checking nextRow
    {
        // CRITICAL REGION STARTED(INCLUDING WHILE CONDITION) ######################
        // here thread gets it's workspace

        size_t row = th->nextRow;
        th->nextRow++; // prevent other threads threads to access this thread working row

        pthread_mutex_unlock (&th->mutex);
        // CRITICAL REGION END #####################################################

        // now the thread can compute the distance matrix inside it's workspace (row)
        register __m256 x1, y1;

        x1 = _mm256_set1_ps(th->inst->X[row]);
        y1 = _mm256_set1_ps(th->inst->Y[row]);

        if (th->inst->params.roundWeights == 0)
        {
            for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
            {
                register __m256 x2 = _mm256_loadu_ps(&th->inst->X[i]);
                register __m256 y2 = _mm256_loadu_ps(&th->inst->Y[i]);

                register __m256 dist =                 

                // store result in memory
                _mm256_storeu_ps(&th->inst->edgeCost.data[row * th->inst->edgeCost.rowElems + i], dist);
            }
            // set last elements of the row (the ones that prevents a lot of headackes with avx) to infinity(so it never gets picked in greedy algorithms)
            for (size_t i = n; i < th->inst->edgeCost.rowElems; i++)
                th->inst->edgeCost.data[row * th->inst->edgeCost.rowElems + i] = INFINITY;
        }
        else
        {
            register __m256i roundedDist;

            for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
            {
                switch (th->inst->params.edgeWeightType)
                {
                case EUC_2D:
                    if (USE_APPROXIMATED_DISTANCES == 1)
                        roundedDist = euclideanCost2DFastApproxRounded(x1, x2, y1, y2);
                    else
                        roundedDist = euclideanCost2DRounded(x1, y1, x2, y2);
                    break;

                case MAN_2D:
                    roundedDist = manhattanCost2DRounded(x1, y1, x2, y2); break;

                case MAX_2D:
                    roundedDist = maximumCost2DRounded(x1, y1, x2, y2); break;

                case ATT:
                    if (USE_APPROXIMATED_DISTANCES == 1)
                        roundedDist = attCost2DFastApproxRounded(x1, y1, x2, y2);
                    else
                        roundedDist = attCost2DRounded(x1, y1, x2, y2);
                    break;

                default:
                    LOG(LOG_LVL_CRITICAL, "Thread is trying to compute distance with unsupported edge weight type. Check the code");
                    return NULL;
                }

                // store result in memory
                _mm256_store_ps((float*)&th->inst->edgeCost.roundedMat[row * th->inst->edgeCost.rowElems + i], (__m256)roundedDist);
            }
            // set last elements of the row (the ones that prevents a lot of headackes with avx) to infinity(so it never gets picked in greedy algorithms)
            for (size_t i = n; i < th->inst->edgeCost.rowElems; i++)
                th->inst->edgeCost.roundedMat[row * th->inst->edgeCost.rowElems + i] = INT_MAX;
        }
    }

    // unlock mutex that was locked in the while condition
    pthread_mutex_unlock(&th->mutex);

    end =  clock();

    double *cpuTimeUsed = malloc(sizeof(double));
    *cpuTimeUsed = ((double)(end-start)) / CLOCKS_PER_SEC;

    pthread_exit(cpuTimeUsed);
}
    