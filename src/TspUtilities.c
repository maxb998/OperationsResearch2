#include "TspUtilities.h"
#include "EdgeCostFunctions.h"

#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h> // used for logger va_list
#include <getopt.h> // args parsing by POSIX
#include <math.h>

#include <immintrin.h>


const char * logLevelString [] = {
	"\033[1;31mERR \033[0m", // 0
	"\033[1;35mCRIT\033[0m", // 1
	"\033[1;33mWARN\033[0m", // 2
	"\033[1;36mNOTI\033[0m", // 3
	"\033[1;34mLOG \033[0m", // 4
	"\033[1;32mDEBG\033[0m", // 5
	"\033[1;90mALL \033[0m"	// 6
};


// Returns the number of processors of the machine
static inline int nProcessors();

Instance newInstance ()
{
    Instance d = { .nNodes = 0, .X = NULL, .Y = NULL, .edgeCostMat = NULL, .params = { .edgeWeightType = 0, .randomSeed = -1, .roundWeightsFlag = 0, .nThreads = nProcessors() } };
    memset(d.params.inputFile, 0, sizeof(d.params.inputFile));
    memset(d.params.name, 0, sizeof(d.params.name));

    return d;
}

Solution newSolution (Instance *inst)
{
    Solution s = { .bestCost = INFINITY, .execTime = 0., .instance = inst };

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
    s.bestCost = sol->bestCost;

    for (size_t i = 0; i < (s.instance->nNodes + AVX_VEC_SIZE * 2); i++)
        s.X[i] = sol->X[i];
    
    for (size_t i = 0; i < s.instance->nNodes + 1; i++)
        s.indexPath[i] = sol->indexPath[i];
    
    return s;
}


int LOG (enum logLevel lvl, char * line, ...)
{
    // check log level
    if (lvl > LOG_LEVEL) return 0;

    // print log level
    printf("[%s] ", logLevelString[lvl]);

    // print passed message and values
    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);
    
    // add new line at the end
    if (line[strlen(line)-1] != '\n')
        printf("\n");

    return 0;
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

void parseArgs (Instance *inst, int argc, char *argv[])
{
    static int roundWeightsFlag = 0, gnuplotFlag = 1;
    static struct option options[] = 
        {
            {"seed", required_argument, 0, 's'},
            {"s", required_argument, 0, 's'},
            {"file", required_argument, 0, 'f'},
            {"f", required_argument, 0, 'f'},
            {"threads", required_argument, 0, 't'},
            {"t", required_argument, 0, 't'},
            {"roundweigths", no_argument, &roundWeightsFlag, 1},
            {"noplot", no_argument, &gnuplotFlag, 0},
            {"output", required_argument, NULL, 'o'},
            {"o", required_argument, NULL, 'o'},
            {"", required_argument, NULL, 'o'},
            {0, 0, 0, 0}
        };
    // set flags in d
    inst->params.roundWeightsFlag = roundWeightsFlag;


    int option_index = 0;
    int opt;

    while ((opt = getopt_long(argc, argv, "s:f:t:o:", options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 's':
            inst->params.randomSeed = (int)strtol(optarg, NULL, 10);
            break;

        case 'f':
            if (access(optarg, R_OK) != 0)
                throwError(inst, NULL, "ERROR: File \"%s\" not found\n", optarg);

            strncpy(inst->params.inputFile, optarg, strlen(optarg)+1);
            break;
        
        case 't':
            inst->params.nThreads = strtoul(optarg, NULL, 10);
            break;
        
        default:
            abort();
        }
    }

    LOG(LOG_LVL_NOTICE, "Received arguments:");
    LOG(LOG_LVL_NOTICE,"    Random Seed  = %d", inst->params.randomSeed);
    LOG(LOG_LVL_NOTICE,"    Filename     = %s", inst->params.inputFile);
    LOG(LOG_LVL_NOTICE,"    Thread Count = %d", inst->params.nThreads);
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

void basicSolutionCheck(Solution *sol)
{
    Instance *inst = sol->instance;

    char * coveredNodes = calloc(sol->instance->nNodes, sizeof(char));

    // First and last node must be equal (the circuit is closed)
    if (sol->indexPath[0] != sol->indexPath[sol->instance->nNodes]) 
        throwError(sol->instance, sol, "BasicSolutionCheck: first and last node in sol.indexPath should coincide, but they do not");

    LOG(LOG_LVL_DEBUG, "BasicSolutionCheck: first and last node in sol.indexPath coincide.");

    // Populate uncoveredNodes array, here we check if a node is repeated along the path
    for (int i = 0; i < inst->nNodes; i++)
    {
        int currentNode = sol->indexPath[i];

        if(coveredNodes[currentNode] == 1)
            throwError(inst, sol, "BasicSolutionCheck: node %d repeated in the solution. Loop iteration %d", currentNode, i);
        else
            coveredNodes[currentNode] = 1;
    }
    LOG(LOG_LVL_EVERYTHING, "BasicSolutionCheck: all nodes in the path are unique.");

    // Check that all the nodes are covered in the path
    for (int i = 0; i < inst->nNodes; i++)
    {
        if(coveredNodes[i] == 0)
            throwError(inst, sol, "BasicSolutionCheck: node %d is not in the path", i);
    }
    free(coveredNodes);
    LOG(LOG_LVL_EVERYTHING, "BasicSolutionCheck: all the nodes are present in the path");


    LOG(LOG_LVL_DEBUG, "BasicSolutionCheck: solution described by sol.indexPath is feasible");

    costCheck(sol);
}

void fullSolutionCheck(Solution *sol)
{
    Instance *inst = sol->instance;

    // also check for sol.X and sol.Y
    if (sol->X[0] != sol->X[inst->nNodes])
        throwError(inst, sol, "FullSolutionCheck: first and last node in sol.X should coincide, but they do not");
    if (sol->Y[0] != sol->Y[inst->nNodes])
        throwError(inst, sol, "FullSolutionCheck: first and last node in sol.Y should coincide, but they do not");

    LOG(LOG_LVL_DEBUG, "FullSolutionCheck: first and last node in sol.X and sol.Y coincide.");

    // Check sol.X and sol.Y if they correspond correctly to inst.X and inst.Y given the indexes in sol.indexPath
    for (size_t i = 0; i < inst->nNodes; i++)
    {
        if (sol->X[i] != inst->X[sol->indexPath[i]])
            throwError(inst, sol, "FullSolutionCheck: sol.X[sol.indexPath[%ld]] = %.3e and does not correspond with inst.X[%ld] = %.3e", i, inst->X[sol->indexPath[i]], i, sol->X[i]);
        if (sol->Y[i] != inst->Y[sol->indexPath[i]])
            throwError(inst, sol, "FullSolutionCheck: sol.Y[sol.indexPath[%ld]] = %.3e and does not correspond with inst.Y[%ld] = %.3e", i, inst->Y[sol->indexPath[i]], i, sol->Y[i]);
    }

    LOG(LOG_LVL_DEBUG, "FullSolutionCheck: solution is coherent and feasible");
}

void costCheck(Solution *sol)
{
    double recomputedCost = computeSolutionCost(sol);

    if (recomputedCost != sol->bestCost)
        throwError(sol->instance, sol, "costChek: Error in the computation of the pathCost. recomputedCost: %lf Cost in Solution: %lf", recomputedCost, sol->bestCost);
    
    LOG(LOG_LVL_DEBUG, "costCheck: Cost of solution is correct");
}

int checkSolutionIntegrity(Solution *sol)
{
    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    for (size_t i = 0; i < n; i++)
    {
        int index = sol->indexPath[i];
        if (index < 0 || index > n)
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: sol.indexPath[%lu] = %d which is not within the limits", i, index);
            return 1;
        }
        if (sol->X[i] != inst->X[index] || sol->Y[i] != inst->Y[index])
        {
            LOG(LOG_LVL_CRITICAL, "checkSolutionIntegrity: Mismatch at index %lu in solution", i);
            return 1;
        }
    }

    // everything checks out
    return 0;
}

void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize, const int printIndex)
{
    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"%s\"\n", sol->instance->params.name);
    fprintf(gnuplotPipe, "set terminal qt size %s\n", plotPixelSize);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pt 7 pointsize %d\n", pointColor, pointSize);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pointsize 0\n", tourPointColor);//, pointSize);


    // assign number to points
    if (printIndex)
        for (size_t i = 0; i < sol->instance->nNodes; i++)
            fprintf(gnuplotPipe, "set label \"%ld\" at %f,%f\n", i, sol->instance->X[i], sol->instance->Y[i]);

    // populating the plot
    
    fprintf(gnuplotPipe, "plot '-' with point linestyle 1, '-' with linespoint linestyle 2\n");

    // first plot only the points
    for (size_t i = 0; i < sol->instance->nNodes; i++)
        fprintf(gnuplotPipe, "%f %f\n", sol->instance->X[i], sol->instance->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    // second print the tour
    for (int i = 0; i <= sol->instance->nNodes; i++)
    {
        fprintf(gnuplotPipe, "%f %f\n", sol->X[i], sol->Y[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    // force write on stream
    fflush(gnuplotPipe);

    // close stream
    pclose(gnuplotPipe);
}



float computeSolutionCostVectorizedFloat(Solution *sol)
{
    register __m256 costVec = _mm256_setzero_ps();
    size_t i = 0;
    while (i < sol->instance->nNodes - AVX_VEC_SIZE)
    {
        register __m256 x1, x2, y1, y2;
        x1 = _mm256_loadu_ps(&sol->X[i]), y1 = _mm256_loadu_ps(&sol->Y[i]);
        x2 = _mm256_loadu_ps(&sol->X[i+1]), y2 = _mm256_loadu_ps(&sol->Y[i+1]);

        costVec = _mm256_add_ps(costVec, squaredEdgeCost_VEC(x1, y1, x2, y2, sol->instance->params.edgeWeightType));

        i += AVX_VEC_SIZE;
    }

    // here we do the last iteration. since all "extra" elements at the end of X and Y are NaNs they can cause issues on sum. so we use a maskload for the last vector
    register __m256 x1, x2, y1, y2;
    x1 = _mm256_loadu_ps(&sol->X[i]); y1 = _mm256_loadu_ps(&sol->Y[i]);
    x2 = _mm256_loadu_ps(&sol->X[i+1]); y2 = _mm256_loadu_ps(&sol->Y[i+1]);

    register __m256 lastDist = squaredEdgeCost_VEC(x1, y1, x2, y2, sol->instance->params.edgeWeightType);

    // now we sustitute the infinity and NaN in the vector with zeroes
    register __m256 infinityVec = _mm256_set1_ps(INFINITY);
    register __m256 mask = _mm256_cmp_ps(lastDist, infinityVec, _CMP_LT_OQ);
    costVec = _mm256_add_ps(costVec, _mm256_blendv_ps(mask, squaredEdgeCost_VEC(x1, y1, x2, y2, sol->instance->params.edgeWeightType), mask));


    float costVecStore[8];
    _mm256_storeu_ps(costVecStore, costVec);
    
    double totalCost = 0.0;
    for (size_t i = 0; i < AVX_VEC_SIZE; i++)
        totalCost += costVecStore[i];
    
    return totalCost;
}

double computeSolutionCostVectorizedDouble(Solution *sol)
{
    register __m256d costVec = _mm256_setzero_pd();
    size_t i = 0;
    while (i < sol->instance->nNodes - AVX_VEC_SIZE)
    {
        register __m256 x1, x2, y1, y2;
        x1 = _mm256_loadu_ps(&sol->X[i]), y1 = _mm256_loadu_ps(&sol->Y[i]);
        x2 = _mm256_loadu_ps(&sol->X[i+1]), y2 = _mm256_loadu_ps(&sol->Y[i+1]);

        register __m256 costVecFloat = squaredEdgeCost_VEC(x1, y1, x2, y2, sol->instance->params.edgeWeightType);

        // convert vector of floats into 2 vectors of doubles and add them to the total cost
        // first half of the vector
        register __m128 partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 0);
        register __m256d partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
        costVec = _mm256_add_pd(costVec, partOfCostVecDouble);
        // second half of the vector
        partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 1);
        partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
        costVec = _mm256_add_pd(costVec, partOfCostVecDouble);

        i += AVX_VEC_SIZE;
    }

    // here we do the last iteration. since all "extra" elements at the end of X and Y are NaNs they can cause issues on sum. so we use a maskload for the last vector
    register __m256 x1, x2, y1, y2;
    x1 = _mm256_loadu_ps(&sol->X[i]); y1 = _mm256_loadu_ps(&sol->Y[i]);
    x2 = _mm256_loadu_ps(&sol->X[i+1]); y2 = _mm256_loadu_ps(&sol->Y[i+1]);

    register __m256 lastDist = squaredEdgeCost_VEC(x1, y1, x2, y2, sol->instance->params.edgeWeightType);

    // now we sustitute the NaNs in the vector with zeroes
    register __m256 infinityVec = _mm256_set1_ps(INFINITY);
    register __m256 mask = _mm256_cmp_ps(lastDist, infinityVec, _CMP_LT_OQ);
    register __m256 costVecFloat = _mm256_blendv_ps(mask, squaredEdgeCost_VEC(x1, y1, x2, y2, sol->instance->params.edgeWeightType), mask);

    // convert vector of floats into 2 vectors of doubles and add them to the total cost
    // first half of the vector
    register __m128 partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 0);
    register __m256d partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
    costVec = _mm256_add_pd(costVec, partOfCostVecDouble);
    // second half of the vector
    partOfCostVecFloat = _mm256_extractf128_ps(costVecFloat, 1);
    partOfCostVecDouble = _mm256_cvtps_pd(partOfCostVecFloat);
    costVec = _mm256_add_pd(costVec, partOfCostVecDouble);

    double vecStore[4];
    _mm256_storeu_pd(vecStore, costVec);
    
    double totalCost = 0.0;
    for (size_t i = 0; i < 4; i++)
        totalCost += vecStore[i];
    
    return totalCost;
}

double computeSolutionCost(Solution *sol)
{
    double cost = 0.0;
    for (size_t i = 0; i < sol->instance->nNodes; i++)
        cost += (double)squaredEdgeCost(sol->X[i], sol->Y[i], sol->X[i+1], sol->Y[i+1], sol->instance->params.edgeWeightType);
    
    return cost;
}