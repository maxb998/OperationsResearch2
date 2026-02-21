#include "Tsp.h"
#include "ArgParser.h"

#include <stdio.h>


const char *modeDoc = "Specify the type of solver to use";
const char *modeNames[] = {"nn", "em", "tabu", "vns", "anneal", "gene", "bnd", "bc", "hf", "lb"};
const char *modeSubDocs[] = {
    "Nearest Neighbor",
    "Extra Mileage",
    "Tabu Search",
    "Variable Neighborhood Search",
    "Simulated Annealing",
    "Genetic Algorithm",
    "Benders method (CPLEX)",
    "Branch & Cut (CPLEX)",
    "Hard Fixing (CPLEX)",
    "Local Branching (CPLEX)"
};

const char *metaInitDoc = "Specify which initialization method to use before running a metaheuristic (DEFAULT=nn)";
#define  META_INIT_COUNT 2

const char *cplexInitDoc = "Specify which initialization method to warm start CPLEX";
#define CPLEX_INIT_COUNT 6

const char *graspDoc = "Specify to Grasp mode (DEFAULT=almostbest)";
const char *graspNames[] = { "none", "almostbest", "random" };
const char *graspSubDocs[] = {
    "Do not use GRASP",
    "The Use of grasp will be limited to selecting another good choice with default probability value",
    "At every iteration have a completely random choice with default probability"
};

const char *logLevelDoc = "Specify the log level (DEFAULT=log)";
const char *logLevelNames[] = { "error", "critical", "warning", "notice", "info", "debug", "trace" };
const char *logLevelSubDocs[] = {
    "Show only error messages",
    "Show critical messages and all above",
    "Show warning and all above",
    "Show notice messages and all above",
    "Show info messages and all above",
    "Show debug messages and all above",
    "Show all messages"
};

const char *compTypeDoc = "Specify the way cost computations are performed";
const char *compTypeNames[] = {"base", "matrix", "avx"};
const char *compTypeSubDocs[] = {
    "Standard way of perfoming cost computations one at a time",
    "Use a matrix to precompute all the costs in the beginning only",
    "Use a avx instructions to perform cost computations where possible"
};


// perform 2 exp by bitshifting
static void intExp2(char *arg, void *dataPtr, void *funcData);


void argParse(Parameters *p, int argc, char *argv[])
{
    ArgGroup groups[] = {
        {.priority=0, .description="Required options"},
        {.priority=1, .description="Heuristics options"},
        {.priority=2, .description="Metaheuristics options"},
        {.priority=3, .description="CPLEX options"},
        {.priority=4, .description="Post-heuristic 2/3-opt manual activation"},
        {.priority=5, .description="Output options"},
        {.priority=6, .description="General options"},
        {.priority=-1}
    };

    ArgOption optionNames[] = {
        {.key='f', .name="file", .required=true, .dtype=DTYPE_STRING, .group=&groups[0], .doc="Location of the .tsp file containing the instance to use", .dataPtr=&p->inputFile},
        {.key='t', .name="tlim", .required=true, .dtype=DTYPE_DOUBLE, .group=&groups[0], .doc="Specify time limit for the execution", .dataPtr=&p->tlim},
        {.key='m', .name="mode", .required=true, .dtype=DTYPE_STRING, .group=&groups[0], .doc=modeDoc, .subNames=modeNames, .subDoc=modeSubDocs, .subCount=SIZE_CONST_ARR(modeNames), .dataPtr=&p->mode, .func=intExp2},

        {.key=0, .name="graspType", .dtype=DTYPE_STRING, .group=&groups[1], .doc=graspDoc, .subNames=graspNames, .subDoc=graspSubDocs, .subCount=SIZE_CONST_ARR(graspNames), .dataPtr=&p->graspType},
        {.key=0, .name="graspChance", .dtype=DTYPE_DOUBLE, .group=&groups[1], .doc="Chance of a grasp event to trigger", .dataPtr=&p->graspChance, .func=NULL},
        {.key=0, .name="nnTryall", .dtype=DTYPE_NONE, .group=&groups[1], .doc="Makes Nearest Neighbor start from each node instead of chosing a random one each time", .dataPtr=&p->nnFirstNodeOption},
        {.key=0, .name="emFarthest", .dtype=DTYPE_NONE, .group=&groups[1], .doc="Makes Extra Mileage initialization the farthest nodes each time instead of a random one each time", .dataPtr=&p->emInitOption},

        {.key=0, .name="metaInit", .dtype=DTYPE_STRING, .group=&groups[2], .doc=metaInitDoc, .subNames=modeNames, .subNames=modeNames, .subDoc=modeSubDocs, .subCount=META_INIT_COUNT, .dataPtr=&p->metaheurInitMode, .func=intExp2},
        {.key=0, .name="metaRestartThreshold", .dtype=DTYPE_UINT, .group=&groups[2], .doc="Specify the threshold for non-improving iterations of vns or tabu berfore restarting from best solution", .dataPtr=&p->metaRestartThreshold},
        {.key=0, .name="tabuTenureSize", .dtype=DTYPE_UINT, .group=&groups[2], .doc="Specify how big the tenure should be in Tabu Search", .dataPtr=&p->tabuTenureSize},
        {.key=0, .name="vnsKickSize", .dtype=DTYPE_UINT, .group=&groups[2], .count=2, .doc="Specify the magnitude of the \"kick\" that randomizes the solution in vns. Eg: --vnsKickSize 2,6", .dataPtr=&p->vnsKickSize},
        {.key=0, .name="geneticParams", .dtype=DTYPE_UINT, .group=&groups[2], .count=4, .doc="Specify the sizes of Population, Crossover, Mutation and Reintroduction in that order in the genetic algorithm. Eg: --geneticParams 50,25,25,5", .dataPtr=&p->geneticParams},
        {.key=0, .name="annealTemperature", .dtype=DTYPE_DOUBLE, .group=&groups[2], .doc="Specify temperature exponent for Simulated Annealing procedure (t = 10^exp)", .dataPtr=&p->annealingTemperature},

        {.key=0, .name="cplexInit", .dtype=DTYPE_STRING, .group=&groups[3], .doc=cplexInitDoc, .subNames=modeNames, .subNames=modeNames, .subDoc=modeSubDocs, .subCount=CPLEX_INIT_COUNT, .dataPtr=&p->matheurInitMode, .func=intExp2},
        {.key=0, .name="cplexDisablePatching", .dtype=DTYPE_NONE, .group=&groups[3], .doc="Disable merging of subtours during benders and branch and cut to build feasible solutions", .dataPtr=&p->cplexPatching, .func=NULL},
        {.key=0, .name="cplexEnableWarmStart", .dtype=DTYPE_NONE, .group=&groups[3], .doc="Enables warm start in CPLEX method using heuristics and metaheuristics beforehand", .dataPtr=&p->cplexWarmStart},
        {.key=0, .name="cplexDisableSolPosting", .dtype=DTYPE_NONE, .group=&groups[3], .doc="Disable cplex posting of solutions during the branch and cut method", .dataPtr=&p->cplexSolPosting, .func=NULL},
        {.key=0, .name="cplexDisableUsercuts", .dtype=DTYPE_NONE, .group=&groups[3], .doc="Disable concorde's functions", .dataPtr=&p->cplexUsercuts, .func=NULL},

        {.key='2', .name="2opt", .dtype=DTYPE_NONE, .group=&groups[4], .doc="Specify to use 2-opt at the end of the selected heuristic", .dataPtr=&p->use2Opt},
        {.key='3', .name="3opt", .dtype=DTYPE_NONE, .group=&groups[4], .doc="Specify to use 3-opt at the end of the selected heuristic", .dataPtr=&p->use3Opt},

        {.key='p', .name="plot", .dtype=DTYPE_NONE, .group=&groups[5], .doc="Specify this if yout want to plot final result", .dataPtr=&p->showPlot},
        {.key='s', .name="save", .dtype=DTYPE_STRING, .group=&groups[5], .doc="Specify this if yout want to save final result in run", .dataPtr=&p->saveSolution},

        {.key=0, .name="seed", .dtype=DTYPE_UINT, .group=&groups[6], .doc="Random Seed [0,MAX_INT32] to use as random seed for the current run. If -1 seed will be random", .dataPtr=&p->randomSeed},
        {.key='j', .name="threads", .dtype=DTYPE_UINT, .group=&groups[6], .doc="Maximum number of threads to use. If not specified gets maximum automatically", .dataPtr=&p->nThreads},
        {.key='r', .name="roundcosts", .dtype=DTYPE_NONE, .group=&groups[6], .doc="Specify this if yout want to use rounded version of edge cost", .dataPtr=&p->roundWeights},
        {.key='l', .name="loglvl", .dtype=DTYPE_STRING, .group=&groups[6], .doc=logLevelDoc, .subNames=logLevelNames, .subDoc=logLevelSubDocs, .subCount=SIZE_CONST_ARR(logLevelNames), .dataPtr=&p->logLevel},
        {.key='c', .name="computationtype", .dtype=DTYPE_STRING, .group=&groups[6], .doc=compTypeDoc, .subNames=compTypeNames, .subDoc=compTypeSubDocs, .subCount=SIZE_CONST_ARR(compTypeNames), .dataPtr=&p->compType},
        {.key=-1} // end
    };

    parser(optionNames, groups, argc, argv);
}

static void intExp2(char *arg, void *dataPtr, void *funcData)
{
    *(int*)dataPtr = 1 << *(int*)dataPtr;
}
