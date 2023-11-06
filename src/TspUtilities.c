#include "Tsp.h"

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h> // used for logger va_list
//#include <math.h>
#include <stdint.h>

const char * logLevelString [] = {
	"\033[1;31mERR \033[0m", // 0
	"\033[1;35mCRIT\033[0m", // 1
	"\033[1;33mWARN\033[0m", // 2
	"\033[1;36mNOTI\033[0m", // 3
	"\033[1;34mLOG \033[0m", // 4
	"\033[1;32mDEBG\033[0m", // 5
	"\033[1;90mALL \033[0m"	// 6
};

static int nProcessors();

Instance newInstance ()
{
    Instance d = {
        .nNodes = 0,
        .X = NULL, .Y = NULL,
        .edgeCostMat = NULL,
        .params = {
            .inputFile = { 0 },
            .mode=MODE_NONE,
            .tlim=-1.,
            
            .graspType=GRASP_NONE,
            .nnFirstNodeOption = NN_FIRST_RANDOM,
            .emInitOption = EM_INIT_RANDOM,

            .metaheurInitMode = MODE_NN,
            .metaRestartThreshold = 1000,
            .tabuTenureSize = -1,
            .vnsKickSize = { .Max=20, .Min=5},
            .geneticParams = { .populationSize=50, .crossoverAmount=25, .mutationAmount=25 },
            
            .matheurInitMode = MODE_NN,
            .hardFixPolicy = HARDFIX_POLICY_RANDOM,

            .use2Opt=false,

            .randomSeed = -1,
            .nThreads = nProcessors(),
            .roundWeights = false,
            .showPlot = false,
            .saveSolution = false,
            .logLevel=LOG_LVL_LOG,

            .edgeWeightType  = -1,
            .name = { 0 },
            .graspChance = 0.1
        }
    };

    return d;
}

static inline int nProcessors()
{
    FILE *commandPipe;
    char *command = "nproc";
    char temp[10];
    commandPipe = (FILE*)popen(command, "r");
    (void)((int*)fgets(temp, sizeof(temp), commandPipe)+1); // evil warning anihilation
    //fgets(temp, sizeof(temp), commandPipe) != NULL ? printf(" \b") : printf(" \b");
    pclose(commandPipe);
    int numProcessors = atoi(temp);
    return numProcessors;
}

Solution newSolution (Instance *inst)
{
    Solution s = {
        .cost = -1LL,
        .execTime = 0.,
        .instance = inst 
    };

    s.indexPath = malloc((inst->nNodes + AVX_VEC_SIZE) * sizeof(int)); // if changing the size of the malloc also change cloneSolution()
    if (!s.indexPath)
        throwError("newSolution: Failed to allocate memory");
    return s;
}

void destroyInstance (Instance *inst)
{
    // memory has been allocated to X, no need to free Y
    free(inst->X);
    free(inst->edgeCostMat);

    inst->X = NULL;
    inst->Y = NULL;
    inst->edgeCostMat = NULL;
}

void destroySolution (Solution *sol)
{
    free(sol->indexPath);
    sol->indexPath = NULL;
}

void cloneSolution(Solution *src, Solution *dst)
{
    dst->cost = src->cost;
    dst->instance = src->instance;
    dst->execTime = src->execTime;
    
    for (int i = 0; i < dst->instance->nNodes + AVX_VEC_SIZE; i++)
        dst->indexPath[i] = src->indexPath[i];
}

static enum LogLevel LOG_LEVEL = LOG_LVL_LOG;
void setLogLevel(enum LogLevel lvl)
{
    LOG_LEVEL = lvl;
}

void LOG (enum LogLevel lvl, char * line, ...)
{
    // check log level
    if (lvl > LOG_LEVEL) return;

    // print log level
    printf("  [%s] ", logLevelString[lvl]);

    // print passed message and values
    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);
    
    // add new line at the end
    if (line[strlen(line)-1] != '\n')
        printf("\n");
}

void throwError (char * line, ...)
{
    printf("[%s] ", logLevelString[0]);

    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);

    printf("\n");

    exit(EXIT_FAILURE);
}

bool checkSolution(Solution *sol)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    char * coveredNodes = calloc(sol->instance->nNodes, sizeof(char));

    // Populate uncoveredNodes array, here we check if a node is repeated along the path
    for (int i = 0; i < n; i++)
    {
        int currentNode = sol->indexPath[i];

        if ((currentNode < 0) || (currentNode > n))
        { LOG(LOG_LVL_CRITICAL, "SolutionCheck: iterpath[%d]=%d is not both >= 0 and < nNodes", i, currentNode); return false; }

        if(coveredNodes[currentNode] == 1)
        { LOG(LOG_LVL_CRITICAL, "SolutionCheck: node %d repeated in the solution. Loop iteration %d", currentNode, i); return false; }
        else
            coveredNodes[currentNode] = 1;
    }

    // Check that all the nodes are covered in the path
    for (int i = 0; i < inst->nNodes; i++)
        if(coveredNodes[i] == 0)
        { LOG(LOG_LVL_CRITICAL, "SolutionCheck: node %d is not in the path", i); return false; }
    free(coveredNodes);

    __uint128_t recomputedCost = computeSolutionCost(sol);

    if (recomputedCost != sol->cost)
    { LOG(LOG_LVL_CRITICAL, "SolutionCheck: Error in the computation of the pathCost. Recomputed Cost: %f Cost in Solution: %f", cvtCost2Double(recomputedCost), cvtCost2Double(sol->cost)); return false; }

    return true;
}

__uint128_t computeSolutionCost(Solution *sol)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    __uint128_t cost = 0;

    #if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
        for (int i = 0; i < n - 1; i++)
            cost += cvtFloat2Cost(computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i+1]], inst->Y[sol->indexPath[i+1]], inst));
        
        cost += cvtFloat2Cost(computeEdgeCost(inst->X[sol->indexPath[n-1]], inst->Y[sol->indexPath[n-1]], inst->X[sol->indexPath[0]], inst->Y[sol->indexPath[0]], inst));
    #elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        for (int i = 0; i < n - 1; i++)
            cost += cvtFloat2Cost(inst->edgeCostMat[(size_t)sol->indexPath[i] * (size_t)n + (size_t)sol->indexPath[i+1]]);
        
        cost += cvtFloat2Cost(inst->edgeCostMat[(size_t)sol->indexPath[n-1] * (size_t)n + (size_t)sol->indexPath[0]]);
    #endif

    return cost;
}

static void quicksort_internal(float *arr, int low, int high)
{
    if ((high < SIZE_MAX) && (low < high))
    {
        int pivot = (int)(((long)rand() + 1) * (long)(high - low) / (long)RAND_MAX) + low;
        swapElems(arr[high], arr[pivot])

        int i = (low - 1);

        for (int j = low; j <= high - 1; j++)
        {
            if (arr[j] <= arr[high])
            {
                i++;
                swapElems(arr[i], arr[j])
            }
        }

        swapElems(arr[i + 1], arr[high])

        quicksort_internal(arr, low, i);
        quicksort_internal(arr, i + 2, high);
    }
}

void sort(float *arr, int n)
{
    quicksort_internal(arr, 0, n-1);
}

static void quickargsort_internal(float *arr, int *indexes, int low, int high)
{
    if ((high < SIZE_MAX) && (low < high))
    {
        int pivot = (int)(((long)rand() + 1) * (long)(high - low) / (long)RAND_MAX) + low;
        swapElems(indexes[high], indexes[pivot])

        int i = (low - 1);
        for (int j = low; j <= high - 1; j++)
        {
            if (arr[indexes[j]] <= arr[indexes[high]])
            {
                i++;
                swapElems(indexes[i], indexes[j])
            }
        }

        swapElems(indexes[i + 1], indexes[high])

        quickargsort_internal(arr, indexes, low, i);
        quickargsort_internal(arr, indexes, i + 2, high);
    }
}

void argsort(float *arr, int *indexes, int n)
{
    for (int i = 0; i < n; i++)
        indexes[i] = i;

    quickargsort_internal(arr, indexes, 0, n-1);
}