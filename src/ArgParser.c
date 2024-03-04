#include "Tsp.h"

#include <argp.h>
#include <unistd.h>
#include <string.h>

#define SUBOPT_BLANKSPACE " "

#define SUBOPT_NN "nn"
#define SUBOPT_EM "em"
#define SUBOPT_TABU "tabu"
#define SUBOPT_VNS "vns"
#define SUBOPT_ANNEALING "annealing"
#define SUBOPT_GENETIC "genetic"
#define SUBOPT_BENDERS "benders"
#define SUBOPT_BRANCHCUT "branch-cut"
#define SUBOPT_HARDFIX "hardfix"
#define SUBOPT_LOCALBRANCHING "local-branching"

#define SUBOPT_GRASP_ALMOSTBEST "almostbest"
#define SUBOPT_GRASP_RANDOM "random"

#define SUBOPT_LOG_ERROR "error"
#define SUBOPT_LOG_CRITICAL "critical"
#define SUBOPT_LOG_WARNING "warning"
#define SUBOPT_LOG_NOTICE "notice"
#define SUBOPT_LOG_LOG "log"
#define SUBOPT_LOG_DEBUG "debug"
#define SUBOPT_LOG_ALL "all"

#define DOC_NN SUBOPT_BLANKSPACE SUBOPT_NN "\t\t: Use Nearest Neighbor\n"
#define DOC_EM SUBOPT_BLANKSPACE SUBOPT_EM "\t\t: Use Extra Mileage\n"
#define DOC_TABU SUBOPT_BLANKSPACE SUBOPT_TABU "\t\t: Use Tabu Search\n"
#define DOC_VNS SUBOPT_BLANKSPACE SUBOPT_VNS "\t\t: Use Variable Neighborhood Search\n"
#define DOC_ANNEALING SUBOPT_BLANKSPACE SUBOPT_ANNEALING "\t\t: Use Simulated Annealing\n"
#define DOC_GENETIC SUBOPT_BLANKSPACE SUBOPT_GENETIC "\t\t: Use Genetic Algorithm\n"
#define DOC_BENDERS SUBOPT_BLANKSPACE SUBOPT_BENDERS "\t: Use Benders method\n"
#define DOC_BRANCHCUT SUBOPT_BLANKSPACE SUBOPT_BRANCHCUT "\t: Use cplex generic callback to perform Branch & Cut\n"
#define DOC_HARDFIX SUBOPT_BLANKSPACE SUBOPT_HARDFIX "\t: Use Hard Fixing Matheuristic\n"
#define DOC_LOCALBRANCHING SUBOPT_BLANKSPACE SUBOPT_LOCALBRANCHING "\t: Use Local Branching Matheuristic\n"

#define ARGP_MODE_DOC "\
Specify the type of solver to use \n"\
DOC_NN DOC_EM DOC_TABU DOC_VNS DOC_ANNEALING DOC_GENETIC DOC_BENDERS DOC_BRANCHCUT DOC_HARDFIX DOC_LOCALBRANCHING

#define HEURISTICS_MODES_COUNT 3
#define METAHEUR_MODES_COUNT 3
#define CPLEX_SOLVERS_COUNT 4

#define MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT + CPLEX_SOLVERS_COUNT)
static const char *modeStrings[] = {
    SUBOPT_NN,
    SUBOPT_EM,
    SUBOPT_TABU,
    SUBOPT_VNS,
    SUBOPT_ANNEALING,
    SUBOPT_GENETIC,
    SUBOPT_BENDERS,
    SUBOPT_BRANCHCUT,
    SUBOPT_HARDFIX,
    SUBOPT_LOCALBRANCHING
};

#define GRASP_DOC "\
Specify to Grasp mode (DEFAULT=random)\n" \
SUBOPT_BLANKSPACE SUBOPT_GRASP_ALMOSTBEST "\t: The Use of grasp will be limited to selecting another good choice with default probability value\n" \
SUBOPT_BLANKSPACE SUBOPT_GRASP_RANDOM "\t: At every iteration have a completely random choice with default probability\n"

static const char *graspStrings[] = { SUBOPT_GRASP_ALMOSTBEST, SUBOPT_GRASP_RANDOM };

#define METAHEURISTICS_INIT_MODE_DOC "\
Specify the heuristic that vns shall use at the start when finding the base solution (DEFAULT=nn)\n" \
DOC_NN DOC_EM
#define METAHEURISTICS_INIT_MODES_COUNT HEURISTICS_MODES_COUNT

#define MATHEUR_INIT_MODE_DOC "\
Specify which Heuristic/Metaheuristic to use as initialization for cplex\n" \
DOC_NN DOC_EM DOC_TABU DOC_VNS DOC_ANNEALING DOC_GENETIC
#define MATHEUR_INIT_MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT)

#define LOG_LEVEL_DOC "\
Specify the log level (DEFAULT=log)\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_ERROR "\t\t: Show only error messages\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_CRITICAL "\t: Show critical messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_WARNING "\t: Show warning and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_NOTICE "\t: Show notice messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_LOG "\t\t: Show log messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_DEBUG "\t\t: Show debug messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_ALL "\t\t: Show all messages\n"

static const char *logLevelStrings[] = { SUBOPT_LOG_ERROR, SUBOPT_LOG_CRITICAL, SUBOPT_LOG_WARNING, SUBOPT_LOG_NOTICE, SUBOPT_LOG_LOG, SUBOPT_LOG_DEBUG, SUBOPT_LOG_ALL };
static const int loglvlsCount = sizeof(logLevelStrings)/sizeof(*logLevelStrings);

enum argpKeys{
    ARGP_FILE='f',
    ARGP_MODE='m',
    ARGP_TLIM='t',

    ARGP_GRASP_MODE=300,
    ARGP_GRASP_CHANCE,
    ARGP_NN_TRYALL,
    ARGP_EM_FARTHEST,
    
    ARGP_META_INIT_MODE,
    ARGP_RESTART_THRESHOLD,
    ARGP_TABU_TENURESIZE,
    ARGP_VNS_KICKSIZE,
    ARGP_GENETIC_PARAMS,
    ARGP_ANNEAL_TEMP,

    ARGP_CPLEX_INIT_MODE,
    ARGP_CPLEX_WARMSTART,
    ARGP_CPLEX_POSTING,
    ARGP_CPLEX_USERCUTS,
    ARGP_HARDFIX_SMALLEST,

    ARGP_2OPT='2',
    ARGP_3OPT='3',

    ARGP_SEED=299,
    ARGP_NTHREADS='j',
    ARGP_ROUND='r',
    ARGP_PLOT='p',
    ARGP_SAVE='s',
    ARGP_LOG_LEVEL='l'
};

error_t argpParser(int key, char *arg, struct argp_state *state);

static void parseModeOption(char *arg, enum Mode *savePtr, const char **optionsSet, const int from, const int to, const char *optionName);
static void parseEnumOption(char *arg, int *savePtr, const char **optionsSet, const int from, const int to, const char *optionName);

static int parseUint(char *arg, char expectedEndChr, const char *paramName);
static void parseUintList(char*arg, const char separator, int *savePtr, int listLenght, const char *paramName);
static double parseDouble(char *arg, const char *paramName);


static void checkEssentials(Instance *inst);


void argParse(Instance * inst, int argc, char *argv[])
{
    static struct argp_option argpOptions[] = {
        { .name="file", .key=ARGP_FILE, .arg="FILENAME", .flags=0, .doc="Location of the .tsp file containing the instance to use\n", .group=1 },
        { .name="mode", .key=ARGP_MODE, .arg="MODE", .flags=0, .doc=ARGP_MODE_DOC, .group=1 },
        { .name="tlim", .key=ARGP_TLIM, .arg="FLOAT", .flags=0, .doc="Specify time limit for the execution\n", .group=1 },

        { .name="graspType", .key=ARGP_GRASP_MODE, .arg="STRING", .flags=0, .doc=GRASP_DOC, .group=2 },
        { .name="graspChance", .key=ARGP_GRASP_CHANCE, .arg="FLOAT", .flags=0, .doc="Probability for a \"gras event\" to happen", .group=2 },
        { .name="nnTryall", .key=ARGP_NN_TRYALL, .arg=NULL, .flags=0, .doc="Specify to make Nearest Neighbor start from each node instead of chosing a random one each time\n", .group=2 },
        { .name="emFarthest", .key=ARGP_EM_FARTHEST, .arg=NULL, .flags=0, .doc="Specify to make Extra Mileage initialization the farthest nodes each time instead of a random one each time\n", .group=2 },

        { .name="metaRestartThreshold", .key=ARGP_RESTART_THRESHOLD, .arg="UINT", .flags=0, .doc="Specify the threshold for non-improving iterations of vns or tabu berfore restarting from best solution\n", .group=3 },
        { .name="metaInit", .key=ARGP_META_INIT_MODE, .arg="STRING", .flags=0, .doc=METAHEURISTICS_INIT_MODE_DOC, .group=3 },
        { .name="tabuTenureSize", .key=ARGP_TABU_TENURESIZE, .arg="UINT", .flags=0, .doc="Specify how big the tenure should be in Tabu Search\n", .group=3 },
        { .name="vnsKickSize", .key=ARGP_VNS_KICKSIZE, .arg="UINT,UINT", .flags=0, .doc="Specify the size of the random \"kick\" that randomizes the solution in vns. Eg: --vnsKickSize 2,6\n", .group=3 },
        { .name="geneticParams", .key=ARGP_GENETIC_PARAMS, .arg="UINT,UINT,UINT", .flags=0, .doc="Specify the sizes of Population, Crossover, Mutation and Reintroduction in that order in the genetic algorithm. Eg: --geneticParams 50,25,25\n", .group=3 },
        { .name="annelTemp", .key=ARGP_ANNEAL_TEMP, .arg="UINT", .flags=0, .doc="Specify temperature for Simulated Annealing procedure\n", .group=3 },

        { .name="cplexInit", .key=ARGP_CPLEX_INIT_MODE, .arg="STRING", .flags=0, .doc=MATHEUR_INIT_MODE_DOC, .group=4 },
        { .name="cplexEnableWarmStart", .key=ARGP_CPLEX_WARMSTART, .arg=NULL, .flags=0, .doc="Enable the ability to find a solution by means of heuristic and metaheuristics and use it to \"warm start\" cplex when using benders of branch and cut methods", .group=4 },
        { .name="cplexDisableSolPosting", .key=ARGP_CPLEX_POSTING, .arg=NULL, .flags=0, .doc="Disable the ability of cplex of posting the best feasible solution found at any point during the branch and cut method", .group=4 },
        { .name="cplexDisableUsercuts", .key=ARGP_CPLEX_USERCUTS, .arg=NULL, .flags=0, .doc="Disable the ability of using concorde's functions to find connected components during cplex relaxation and add cuts that violate such components as usercuts", .group=4 },
        { .name="hardfixSmallest", .key=ARGP_HARDFIX_SMALLEST, .arg=NULL, .flags=0, .doc="Specify to make Hard Fixing fix only edges with smallest cost instead of fixing random edges\n", .group=4 },

        { .name="2opt", .key=ARGP_2OPT, .arg=NULL, .flags=0, .doc="Specify to use 2-opt at the end of the selected heuristic\n", .group=5 },
        { .name="3opt", .key=ARGP_3OPT, .arg=NULL, .flags=0, .doc="Specify to use 3-opt at the end of the selected heuristic\n", .group=5 },

        { .name="seed", .key=ARGP_SEED, .arg="UINT", .flags=0, .doc="Random Seed [0,MAX_INT32] to use as random seed for the current run. If -1 seed will be random\n", .group=6 },
        { .name="threads", .key=ARGP_NTHREADS, .arg="UINT", .flags=0, .doc="Maximum number of threads to use. If not specified gets maximum automatically\n", .group=6 },
        { .name="roundcosts", .key=ARGP_ROUND, .arg=NULL, .flags=0, .doc="Specify this if yout want to use rounded version of edge cost\n", .group=6 },
        { .name="plot", .key=ARGP_PLOT, .arg=NULL, .flags=0, .doc="Specify this if yout want to plot final result\n", .group=6 },
        { .name="save", .key=ARGP_SAVE, .arg=NULL, .flags=0, .doc="Specify this if yout want to save final result in run/\n", .group=6 },
        { .name="loglvl", .key=ARGP_LOG_LEVEL, .arg="STRING", .flags=0, .doc=LOG_LEVEL_DOC, .group=6 },
        { 0 }
    };

    static struct argp argpData = {
        .options = argpOptions,
        .parser=argpParser, 
    };

    argp_parse(&argpData, argc, argv, 0, 0, inst);
}

error_t argpParser(int key, char *arg, struct argp_state *state)
{
    Instance *inst = state->input;

    switch (key)
    {
    case ARGP_FILE: // get the input filename
        if (access(arg, R_OK)) // check if file exists and is accessible
        {
            LOG(LOG_LVL_ERROR, "File \"%s\" cannot be accessed or does not exist", arg);
            return ARGP_ERR_UNKNOWN;
        }
        strncpy(inst->params.inputFile, arg, strlen(arg));
        break;

    case ARGP_MODE:
        parseModeOption(arg, &inst->params.mode, modeStrings, 0, MODES_COUNT, "mode");
        break;

    case ARGP_TLIM:
        inst->params.tlim = parseDouble(arg, "tlim");
        if (inst->params.tlim <= 0)
            throwError("Time limit cannot be zero or negative");
        break;

    case ARGP_GRASP_MODE:
        parseEnumOption(arg, (int*)&inst->params.graspType, graspStrings, 0, 2, "graspType");
        break;
    
    case ARGP_GRASP_CHANCE:
        inst->params.graspChance = parseDouble(arg, "graspChance");
        if ((inst->params.graspChance < 0.) || (inst->params.graspChance > 1.))
            throwError("Grasp Chance must be inside the interval [0,1]");
        break;

    case ARGP_NN_TRYALL:
        inst->params.nnFirstNodeOption = NN_FIRST_TRYALL;
        break;
    
    case ARGP_EM_FARTHEST:
        inst->params.emInitOption = EM_INIT_FARTHEST_POINTS;
        break;
    
    case ARGP_META_INIT_MODE:
        parseEnumOption(arg, (int*)&inst->params.metaheurInitMode, modeStrings, 0, HEURISTICS_MODES_COUNT, "metaInit");
        break;

    case ARGP_RESTART_THRESHOLD:
        inst->params.metaRestartThreshold = parseUint(arg, 0, "metaRestartThreshold");
        break;

    case ARGP_TABU_TENURESIZE:
        inst->params.tabuTenureSize = parseUint(arg, 0, "tabuTenureSize");
        break;
    
    case ARGP_VNS_KICKSIZE:
        parseUintList(arg, ',', (int*)&inst->params.vnsKickSize, 2, "vnsKickSize");

        if ((inst->params.vnsKickSize.Max < 2) || inst->params.vnsKickSize.Min < 2)
            throwError("vnsKickSize components must be both positive greater than two");
        if (inst->params.vnsKickSize.Max < inst->params.vnsKickSize.Min)
            throwError("vnsKickSize must be specified in the current format: MIN-KICK,MAX-KICK. Eg. 5,20   MIN_KICK cannot be greater or equal than MAX_KICK");
        break;
    
    case ARGP_GENETIC_PARAMS:
        parseUintList(arg, ',', (int*)&inst->params.geneticParams, 4, "geneticParams");
        break;
    
    case ARGP_ANNEAL_TEMP:
        inst->params.annealingTemperature = parseUint(arg, 0, "annelTemp");
        break;

    case ARGP_CPLEX_INIT_MODE:
        parseModeOption(arg, &inst->params.matheurInitMode, modeStrings, 0, HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT, "cplexInit");
        break;

    case ARGP_CPLEX_WARMSTART:
        inst->params.cplexWarmStart = true;
        break;
    
    case ARGP_CPLEX_POSTING:
        inst->params.cplexSolPosting = false;
        break;

    case ARGP_CPLEX_USERCUTS:
        inst->params.cplexUsercuts = false;
        break;

    case ARGP_HARDFIX_SMALLEST:
        inst->params.hardFixPolicy = HARDFIX_POLICY_SMALLEST;
        break;

    case ARGP_2OPT:
        inst->params.use2Opt = true;
        break;
    
    case ARGP_3OPT:
        inst->params.use3Opt = true;
        break;

    case ARGP_SEED:
        inst->params.randomSeed = parseUint(arg, 0, "seed");
        break;

    case ARGP_NTHREADS:
        inst->params.nThreads = parseUint(arg, 0, "nThreads");
        break;

    case ARGP_ROUND:
        inst->params.roundWeights = true;
        break;

    case ARGP_PLOT:
        inst->params.showPlot = true;
        break;
    
    case ARGP_SAVE:
        inst->params.saveSolution = true;
        break;
    
    case ARGP_LOG_LEVEL:
        parseEnumOption(arg, (int*)&inst->params.logLevel, logLevelStrings, 0, loglvlsCount, "loglvl");
        setLogLevel(inst->params.logLevel);
        break;
    
    case ARGP_KEY_END:
        // check if necessary flags have been provided
        checkEssentials(inst);
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static void parseModeOption(char *arg, enum Mode *savePtr, const char **optionsSet, const int from, const int to, const char *optionName)
{
    for (int i = from; i < to; i++)
    {
        if (strcmp(arg, optionsSet[i]) == 0)
        {
            *savePtr = (int)powl(2, i);
            return;
        }
    }
    
    throwError("%s: argument not valid", optionName);
}

static void parseEnumOption(char *arg, int *savePtr, const char **optionsSet, const int from, const int to, const char *optionName)
{
    for (int i = from; i < to; i++)
    {
        if (strcmp(arg, optionsSet[i]) == 0)
        {
            *savePtr = i;
            return;
        }
    }
    
    throwError("%s: argument not valid", optionName);
}

static int parseUint(char *arg, char expectedEndChr, const char *paramName)
{
    char *endPtr;
    long cvt = strtol(arg, &endPtr, 10);
    if (cvt < 0)
        throwError("The value specified as %s cannot be negative", paramName);
    if (*endPtr != expectedEndChr)
        throwError("There are extra character after the %s value or formatting is not correct. Check formats with --help", paramName);

    return (int)cvt;
}

static void parseUintList(char*arg, const char separator, int *savePtr, int listLenght, const char *paramName)
{
    char *endPtr = arg, *startPtr = arg;
    for (int i = 0; i < listLenght; i++)
    {
        if (startPtr == NULL)
            throwError("Missing and element for the option %s. Check --help", paramName);

        savePtr[i] = (int)strtol(startPtr, &endPtr, 10);
        if (savePtr[i] < 0)
            throwError("Cannot use negative numbers in %s option", paramName);
        if ((*endPtr != separator) && (*endPtr != 0))
            throwError("Option %s not formatted correctly. Check correct format in --help", paramName);

        startPtr = endPtr + 1;
    }
}

static double parseDouble(char *arg, const char *paramName)
{
    char *endPtr;
    double cvt = strtod(arg, &endPtr);
    if (cvt <= 0)
        throwError("The value specified as %s must be a real number", paramName);
    if (endPtr != &arg[strlen(arg)])
        LOG(LOG_LVL_WARNING, "There are extra character after the %s value", paramName);

    return cvt;
}

static void checkEssentials(Instance *inst)
{
    if (inst->params.mode == MODE_NONE)
        throwError("Missing --mode (-m) option in arguments. A mode must be specified. Check --help or --usage to check what modes are avilable");
    if (inst->params.inputFile[0] == 0)
        throwError("Missing --file (-f) option in arguments. An input file must be specified");
    if (inst->params.tlim == -1.)
    {
        if (inst->params.graspType != GRASP_NONE)
            throwError("If grasp is used a time limit must be provided with option --tlim (-t)");
        if (inst->params.mode >= MODE_VNS)
            throwError("When using an heuristic that is intrinsically using grasp(like %s), a time limit must also be specified", modeStrings[inst->params.mode]);
    }
}


void printInfo(Instance *inst)
{
    Parameters *p = &inst->params;

    bool useNN = (p->mode & MODE_NN) || 
        (p->mode & (MODE_TABU | MODE_VNS | MODE_ANNEALING) & (p->metaheurInitMode & MODE_NN)) || 
        (p->mode & (MODE_BENDERS | MODE_BRANCH_CUT) && p->cplexWarmStart && ((p->matheurInitMode & MODE_NN) | (p->metaheurInitMode & MODE_NN) && (p->matheurInitMode & (MODE_TABU | MODE_VNS | MODE_ANNEALING))));
    bool useEM = (p->mode & MODE_EM) || 
        (p->mode & (MODE_TABU | MODE_VNS | MODE_ANNEALING) & (p->metaheurInitMode & MODE_EM)) || 
        (p->mode & (MODE_BENDERS | MODE_BRANCH_CUT) && p->cplexWarmStart && ((p->matheurInitMode & MODE_EM) | (p->metaheurInitMode & MODE_EM) && (p->matheurInitMode & (MODE_TABU | MODE_VNS | MODE_ANNEALING))));

    printf("SETTINGS:\n");

    // input file
    printf("\t" "Input File/Problem: \"%s\"\n", p->inputFile);
    // mode
    printf("\t" "Current running mode is %s\n", modeStrings[(int)log2l(p->mode)]);

    // time limit
    if (p->tlim != -1)
        printf("\tTime limit of %lf seconds\n", p->tlim);
    else
        printf("\tTime limit is not set\n");

    // grasp
    if (useNN | useNN)
    {
        if (p->graspType == GRASP_NONE)
            printf("\tGrasp is off\n");
        else
            printf("\tGrasp is on, grasp mode id is %s with chance %lf\n", graspStrings[p->graspType], p->graspChance);
    }
    // 2Opt
    if (p->use2Opt || (p->mode & (MODE_TABU | MODE_VNS)))
        printf("\tUsing 2Opt\n");
    
    // metaheuristics modes
    if ((p->mode & (MODE_TABU | MODE_VNS | MODE_ANNEALING)) ||
        ((p->mode & (MODE_BENDERS | MODE_BRANCH_CUT)) && p->cplexWarmStart) ||
        (p->mode & (MODE_HARDFIX | MODE_LOCAL_BRANCHING) & p->matheurInitMode & (MODE_TABU | MODE_VNS | MODE_ANNEALING)))
        printf("\tMetaheuristics initialization set to: %s\n", modeStrings[(int)log2l(p->metaheurInitMode)]);
    // matheuristics modes
    if (((p->mode & (MODE_BENDERS | MODE_BRANCH_CUT)) && p->cplexWarmStart) || 
        (p->mode & (MODE_HARDFIX | MODE_LOCAL_BRANCHING)))
        printf("\tMatheuristics/Cplex initialization set to: %s\n", modeStrings[(int)log2l(p->matheurInitMode)]);

    // nn options
    if (useNN)
        printf("\tNearest Neighbor starting node set to: %s\n", inst->params.nnFirstNodeOption == NN_FIRST_RANDOM ? "random" : "tryall");
    // em options
    if (useEM)
        printf("\tExtra Mileage initialization set to: %s\n", inst->params.emInitOption == EM_INIT_RANDOM ? "random" : "farthest");
    // cplex options
    if (p->mode & (MODE_BENDERS | MODE_BRANCH_CUT | MODE_LOCAL_BRANCHING | MODE_HARDFIX))
        if (p->cplexWarmStart)
            printf("\tCplex will be warm started using an heuristic/metaheuristic solution\n");
    if (p->mode & (MODE_BRANCH_CUT | MODE_LOCAL_BRANCHING | MODE_HARDFIX))
    {
        if (!p->cplexSolPosting)
            printf("\tCplex solution posting is disabled\n");
        if (!p->cplexWarmStart)
            printf("\tCplex usercuts are disabled\n");
    }
    // hardfix options
    if (p->mode == MODE_HARDFIX)
        printf("\tHard Fixing policy set to: %s\n", inst->params.hardFixPolicy == HARDFIX_POLICY_RANDOM ? "random" : "smallest");

    // seed
    if (p->randomSeed != -1)
        printf("\tRandom Seed = %d\n", p->randomSeed);
    // threads
    printf("\tThreads used = %d\n", p->nThreads);
    // roundcosts
    if (p->roundWeights)
        printf("\tEdge Cost is rounded according to tsplib documentation file\n");
    // save
    if (p->saveSolution)
        printf("\tFinal solution of this run will be saved in a .tour file inside OperationsResearch2/runs\n");
    // log level
    printf("\tLog level = %s", logLevelStrings[p->logLevel]);

    printf("\n");
}
