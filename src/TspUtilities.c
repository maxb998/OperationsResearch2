#include "Tsp.h"

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h> // used for logger va_list
#include <math.h>

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
            .use2OptFlag=0,
            .tlim=-1.,

            .nnFirstNodeOption = NN_FIRST_RANDOM,
            .emInitOption = EM_INIT_RANDOM,

            .metaheurInitMode = MODE_NN,
            .matheurInitMode = MODE_NN,
            .warmStartMode = MODE_NN,

            .hardFixPolicy = HARDFIX_POLICY_RANDOM,

            .randomSeed = -1,
            .nThreads = nProcessors(),
            .roundWeightsFlag = 0,
            .showPlotFlag=0,
            .saveFlag=0,
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

    // allocate memory (consider locality). Alse leave place at the end to repeat the first element of the solution(some algoritms benefits from it)
    s.X = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float) + (inst->nNodes + 1) * sizeof(int));
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

Solution cloneSolution(Solution *sol)
{
    Solution s = newSolution(sol->instance);
    s.cost = sol->cost;

    for (size_t i = 0; i < (s.instance->nNodes + AVX_VEC_SIZE * 2); i++)
        s.X[i] = sol->X[i];
    
    for (size_t i = 0; i < s.instance->nNodes + 1; i++)
        s.indexPath[i] = sol->indexPath[i];
    
    return s;
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

void checkSolution(Solution *sol)
{
    Instance *inst = sol->instance;

    char * coveredNodes = calloc(sol->instance->nNodes, sizeof(char));

    // First and last node must be equal (the circuit is closed)
    if (sol->indexPath[0] != sol->indexPath[sol->instance->nNodes]) 
        throwError(sol->instance, sol, "SolutionCheck: first and last node in sol.indexPath should coincide, but they do not");
    // also check for sol.X and sol.Y
    if (sol->X[0] != sol->X[inst->nNodes])
        throwError(inst, sol, "SolutionCheck: first and last node in sol.X should coincide, but they do not");
    if (sol->Y[0] != sol->Y[inst->nNodes])
        throwError(inst, sol, "SolutionCheck: first and last node in sol.Y should coincide, but they do not");

    LOG(LOG_LVL_EVERYTHING, "SolutionCheck: first and last node in sol.indexPath coincide.");

    // Populate uncoveredNodes array, here we check if a node is repeated along the path
    for (int i = 0; i < inst->nNodes; i++)
    {
        int currentNode = sol->indexPath[i];

        if(coveredNodes[currentNode] == 1)
            throwError(inst, sol, "SolutionCheck: node %d repeated in the solution. Loop iteration %d", currentNode, i);
        else
            coveredNodes[currentNode] = 1;
    }
    LOG(LOG_LVL_EVERYTHING, "SolutionCheck: all nodes in the path are unique.");

    // Check that all the nodes are covered in the path
    for (int i = 0; i < inst->nNodes; i++)
        if(coveredNodes[i] == 0)
            throwError(inst, sol, "SolutionCheck: node %d is not in the path", i);
    free(coveredNodes);
    LOG(LOG_LVL_EVERYTHING, "SolutionCheck: all the nodes are present in the path");

    // Check sol.X and sol.Y if they correspond correctly to inst.X and inst.Y given the indexes in sol.indexPath
    for (size_t i = 0; i < inst->nNodes; i++)
    {
        if (sol->X[i] != inst->X[sol->indexPath[i]])
            throwError(inst, sol, "SolutionCheck: sol.X[sol.indexPath[%ld]] = %.3e and does not correspond with inst.X[%ld] = %.3e", i, inst->X[sol->indexPath[i]], i, sol->X[i]);
        if (sol->Y[i] != inst->Y[sol->indexPath[i]])
            throwError(inst, sol, "SolutionCheck: sol.Y[sol.indexPath[%ld]] = %.3e and does not correspond with inst.Y[%ld] = %.3e", i, inst->Y[sol->indexPath[i]], i, sol->Y[i]);
    }

    double recomputedCost = computeSolutionCost(sol);//computeSolutionCostVectorized(sol);

    if (recomputedCost != sol->cost)
        throwError(sol->instance, sol, "SolutionCheck: Error in the computation of the pathCost. Recomputed Cost: %lf Cost in Solution: %lf", recomputedCost, sol->cost);
    
    LOG(LOG_LVL_EVERYTHING, "SolutionCheck: Cost of solution is correct");

    LOG(LOG_LVL_DEBUG, "SolutionCheck: solution is coherent and feasible");
}

double computeSolutionCostVectorized(Solution *sol)
{
    size_t n = sol->instance->nNodes;
    enum EdgeWeightType ewt = sol->instance->params.edgeWeightType ;
    int roundFlag = sol->instance->params.roundWeightsFlag;

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
    
    double totalCost = 0.0;
    for (size_t i = 0; i < 4; i++)
        totalCost += vecStore[i];
    
    return totalCost;
}

double computeSolutionCost(Solution *sol)
{
    enum EdgeWeightType ewt = sol->instance->params.edgeWeightType ;
    int roundFlag = sol->instance->params.roundWeightsFlag;

    double cost = 0.0;
    for (size_t i = 0; i < sol->instance->nNodes; i++)
        cost += (double)computeEdgeCost(sol->X[i], sol->Y[i], sol->X[i+1], sol->Y[i+1], ewt, roundFlag);
    
    return cost;
}