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

    // allocate memory (consider locality). Also leave place at the end to repeat the first element of the solution(some algoritms benefits from it)
    s.X = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float) + (inst->nNodes + AVX_VEC_SIZE) * sizeof(int));
    s.Y = &s.X[inst->nNodes + AVX_VEC_SIZE];
    s.indexPath = (int*)&s.Y[inst->nNodes + AVX_VEC_SIZE];

    return s;
}

void destroyInstance (Instance *inst)
{
    // memory has been allocated to X, no need to free Y
    free(inst->X);
    free(inst->edgeCostMat);

    // reset pointers to NULL
    inst->X = inst->Y = inst->edgeCostMat = NULL;
}

void destroySolution (Solution *sol)
{
    // memory has been allocated to X, no need to free Y or indexPath
    free(sol->X);
    //free(sol->indexPath);

    // reset pointers to NULL
    sol->X = sol->Y = NULL;
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

    for (size_t i = 0; i < ((dst->instance->nNodes + AVX_VEC_SIZE) * 2); i++)
        dst->X[i] = src->X[i];
    
    for (size_t i = 0; i < dst->instance->nNodes + 1; i++)
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

    char * coveredNodes = calloc(sol->instance->nNodes, sizeof(char));

    // First and last node must be equal (the circuit is closed)
    if (sol->indexPath[0] != sol->indexPath[sol->instance->nNodes])
    { LOG(LOG_LVL_CRITICAL, "SolutionCheck: first and last node in sol.indexPath should coincide, but they do not"); return false; }
    // also check for sol.X and sol.Y
    if (sol->X[0] != sol->X[inst->nNodes])
    { LOG(LOG_LVL_CRITICAL, "SolutionCheck: first and last node in sol.X should coincide, but they do not"); return false; }
    if (sol->Y[0] != sol->Y[inst->nNodes])
    { LOG(LOG_LVL_CRITICAL, "SolutionCheck: first and last node in sol.Y should coincide, but they do not"); return false; }

    // Populate uncoveredNodes array, here we check if a node is repeated along the path
    for (int i = 0; i < inst->nNodes; i++)
    {
        int currentNode = sol->indexPath[i];

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

    // Check sol.X and sol.Y if they correspond correctly to inst.X and inst.Y given the indexes in sol.indexPath
    for (size_t i = 0; i < inst->nNodes; i++)
    {
        if (sol->X[i] != inst->X[sol->indexPath[i]])
        { LOG(LOG_LVL_CRITICAL, "SolutionCheck: sol.X[sol.indexPath[%ld]] = %.3e and does not correspond with inst.X[%ld] = %.3e", i, inst->X[sol->indexPath[i]], i, sol->X[i]); return false; }
        if (sol->Y[i] != inst->Y[sol->indexPath[i]])
        { LOG(LOG_LVL_CRITICAL, "SolutionCheck: sol.Y[sol.indexPath[%ld]] = %.3e and does not correspond with inst.Y[%ld] = %.3e", i, inst->Y[sol->indexPath[i]], i, sol->Y[i]); return false; }
    }

    double recomputedCost = computeSolutionCost(sol);//computeSolutionCostVectorized(sol);

    if (recomputedCost != sol->cost)
    { LOG(LOG_LVL_CRITICAL, "SolutionCheck: Error in the computation of the pathCost. Recomputed Cost: %lf Cost in Solution: %lf", recomputedCost, sol->cost); return false; }

    return true;
}

double computeSolutionCostVectorized(Solution *sol)
{
    size_t n = sol->instance->nNodes;
    enum EdgeWeightType ewt = sol->instance->params.edgeWeightType ;
    bool roundFlag = sol->instance->params.roundWeights;

    __m256d costVec = _mm256_setzero_pd();
    size_t i = 0;
    while (i < n) //- AVX_VEC_SIZE)
    {
        __m256 x1, x2, y1, y2;

        x1 = _mm256_loadu_ps(&sol->X[i]);
        y1 = _mm256_loadu_ps(&sol->Y[i]);
        x2 = _mm256_loadu_ps(&sol->X[i+1]);
        y2 = _mm256_loadu_ps(&sol->Y[i+1]);

        if ((i > n - AVX_VEC_SIZE) && (n % AVX_VEC_SIZE != 0))
        {
            int loadMask[AVX_VEC_SIZE] = { 0 };
            for (size_t j = 0; j < n % AVX_VEC_SIZE; j++)
                loadMask[j] = -1;
            __m256i mask = _mm256_loadu_si256((__m256i_u*)loadMask);
            
            x1 = _mm256_maskload_ps(&sol->X[i], mask);
            y1 = _mm256_maskload_ps(&sol->Y[i], mask);
            x2 = _mm256_maskload_ps(&sol->X[i+1], mask);
            y2 = _mm256_maskload_ps(&sol->Y[i+1], mask);
        }

        __m256 costVecFloat = computeEdgeCost_VEC(x1, y1, x2, y2, ewt, roundFlag);

        // convert vector of floats into 2 vectors of doubles and add them to the total cost
        // first half of the vector
        __m128 partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 0);
        __m256d partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
        costVec = _mm256_add_pd(costVec, partOfCostVecDouble);
        // second half of the vector
        partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 1);
        partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
        costVec = _mm256_add_pd(costVec, partOfCostVecDouble);

        i += AVX_VEC_SIZE;
    }

    double vecStore[4];
    _mm256_storeu_pd(vecStore, costVec);
    
    double totalCost = 0;
    for (size_t i = 0; i < 4; i++)
        totalCost += vecStore[i];
    
    return totalCost;
}

double computeSolutionCost(Solution *sol)
{
    enum EdgeWeightType ewt = sol->instance->params.edgeWeightType ;
    bool roundFlag = sol->instance->params.roundWeights;

    double cost = 0;
    for (size_t i = 0; i < sol->instance->nNodes; i++)
        cost += (double)computeEdgeCost(sol->X[i], sol->Y[i], sol->X[i+1], sol->Y[i+1], ewt, roundFlag);
    
    return cost;
}

static void quicksort_internal(float *arr, size_t low, size_t high)
{
    if ((high < SIZE_MAX) && (low < high))
    {
        size_t pivot = (((size_t)rand() + 1) * (high - low) / RAND_MAX) + low;
        {
            register float temp;
            swapElems(arr[high], arr[pivot], temp);
        }

        size_t i = (low - 1);

        for (size_t j = low; j <= high - 1; j++)
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

void sort(float *arr, size_t n)
{
    quicksort_internal(arr, 0, n-1);
}

static void quickargsort_internal(float *arr, int *indexes, size_t low, size_t high)
{
    if ((high < SIZE_MAX) && (low < high))
    {
        size_t pivot = (((size_t)rand() + 1) * (high - low) / RAND_MAX) + low;
        {
            register int tempi;
            swapElems(indexes[high], indexes[pivot], tempi);
        }

        size_t i = (low - 1);
        for (size_t j = low; j <= high - 1; j++)
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

void argsort(float *arr, int *indexes, size_t n)
{
    for (size_t i = 0; i < n; i++)
        indexes[i] = (int)i;

    quickargsort_internal(arr, indexes, 0, n-1);
}