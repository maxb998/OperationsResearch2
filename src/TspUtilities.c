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
            
            .graspType=GRASP_NONE,
            .use2Opt=false,
            .tlim=-1.,

            .nnFirstNodeOption = NN_FIRST_RANDOM,
            .emInitOption = EM_INIT_RANDOM,

            .metaheurInitMode = MODE_NN,
            .matheurInitMode = MODE_NN,
            .warmStartMode = MODE_NN,

            .hardFixPolicy = HARDFIX_POLICY_RANDOM,

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
    fgets(temp, sizeof(temp), commandPipe);
    pclose(commandPipe);
    int numProcessors = atoi(temp);
    return numProcessors;
}

Solution newSolution (Instance *inst)
{
    Solution s = { .cost = INFINITY, .execTime = 0., .instance = inst };
    s.indexPath = malloc((inst->nNodes + AVX_VEC_SIZE) * sizeof(int));
    if (!s.indexPath)
        throwError(inst, NULL, "newSolution: Failed to allocate memory");
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
    if (src->instance->nNodes != dst->instance->nNodes)
    {
        destroySolution(dst);
        throwError(src->instance, src, "cloneSolution: Size of src and dst does not match");
    }
    dst->cost = src->cost;
    dst->instance = src->instance;
    dst->execTime = src->execTime;
    
    for (int i = 0; i < dst->instance->nNodes + 1; i++)
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

void throwError (Instance *inst, Solution *sol, char * line, ...)
{
    printf("[%s] ", logLevelString[0]);

    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);

    printf("\n");

    // free allocated memory
    if (inst)
        destroyInstance(inst);
    if (sol)
    {
        destroyInstance(sol->instance);
        destroySolution(sol);
    }

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

    double recomputedCost = computeSolutionCost(sol);

    if (recomputedCost != sol->cost)
    { LOG(LOG_LVL_CRITICAL, "SolutionCheck: Error in the computation of the pathCost. Recomputed Cost: %lf Cost in Solution: %lf", recomputedCost, sol->cost); return false; }

    return true;
}

double computeSolutionCost(Solution *sol)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;
    enum EdgeWeightType ewt = inst->params.edgeWeightType ;
    bool roundFlag = inst->params.roundWeights;

    double cost = 0;
    for (int i = 0; i < n - 1; i++)
        cost += (double)computeEdgeCost(inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]], inst->X[sol->indexPath[i+1]], inst->Y[sol->indexPath[i+1]], ewt, roundFlag);
    
    cost += (double)computeEdgeCost(inst->X[sol->indexPath[n-1]], inst->Y[sol->indexPath[n-1]], inst->X[sol->indexPath[0]], inst->Y[sol->indexPath[0]], ewt, roundFlag);

    return cost;
}

static void quicksort_internal(float *arr, int low, int high)
{
    if ((high < SIZE_MAX) && (low < high))
    {
        int pivot = (int)(((long)rand() + 1) * (long)(high - low) / (long)RAND_MAX) + low;
        {
            register float temp;
            swapElems(arr[high], arr[pivot], temp);
        }

        int i = (low - 1);

        for (int j = low; j <= high - 1; j++)
        {
            if (arr[j] <= arr[high])
            {
                i++;
                register float tempf;
                swapElems(arr[i], arr[j], tempf);
            }
        }

        register float tempf;
        swapElems(arr[i + 1], arr[high], tempf);

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
        {
            register int tempi;
            swapElems(indexes[high], indexes[pivot], tempi);
        }

        int i = (low - 1);
        for (int j = low; j <= high - 1; j++)
        {
            if (arr[indexes[j]] <= arr[indexes[high]])
            {
                i++;
                register int tempi;
                swapElems(indexes[i], indexes[j], tempi);
            }
        }

        register int tempi;
        swapElems(indexes[i + 1], indexes[high], tempi);

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