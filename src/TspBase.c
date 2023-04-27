#include "TspBase.h"

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

            .randomSeed = -1,
            .nThreads = nProcessors(),
            .roundWeightsFlag = 0,
            .showPlotFlag=0,
            .saveFlag=0,
            .logLevel=LOG_LVL_LOG,

            .edgeWeightType = -1,
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

static enum logLevel LOG_LEVEL = LOG_LVL_LOG;
void setLogLevel(enum logLevel lvl)
{
    LOG_LEVEL = lvl;
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