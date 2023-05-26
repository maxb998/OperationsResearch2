#include "Tsp.h"

#include <argp.h>
#include <unistd.h>
#include <string.h>

#define ARGP_MODE_DOC "\
Specify the type of solver to use: \n\
    nn          : Use Nearest Neighbor\n\
    em          : Use Extra Mileage\n\
    tabu        : Use Tabu Search\n\
    vns         : Use Variable Neighborhood Search\
    annealing   : Use Simulated Annealing\n\
    genetic     : Use Genetic Heuristic\n\
    benders     : Use Benders method\n\
    branch_cut  : Use cplex generic callback to perform Branch & Cut\n\
    hardfix     : Use Hard Fixing Matheuristic\n\
    local-branching : Use Local Branching Matheuristic\n\
"

#define HEURISTICS_MODES_COUNT 2
#define METAHEUR_MODES_COUNT 4
#define EXACT_SOLVERS_COUNT 2
#define MATHEURISTICS_COUNT 2

#define MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT + EXACT_SOLVERS_COUNT + MATHEURISTICS_COUNT)
static const char *modeStrings[] = {
    "nn",
    "em",
    "tabu",
    "vns",
    "annealing",
    "genetic",
    "benders",
    "branch_cut",
    "hardfix",
    "local-branching"
};

#define ARGP_GRASP_DOC "\
Specify to Grasp mode (DEFAULT=random(0.1))\n\
    almostbest         : The Use of grasp will be limited to selecting another good choice with default probability value \n\
    almostbest(<GRASP_CHANCE>) : Same as \"almostbest\" but the selection of the secondary choice happens with probability specified with <GRASP_CHANCE>\n\
    random             : At every iteration have a completely random choice with default probability\n\
    random(<GRASP_CHANCE>)     : Same as \"random\" but the probability of the random choice is specified with <GRASP_CHANCE>\n\
"

static const char *graspStrings[] = {
    "almostbest",
    "random"
};

#define ARGP_NN_DOC "\
Specify the way the first node is selected when using Nearest Neighbor anywhere in the heuristic (DEFAULT=random)\n\
    random              : Select a node at random\n\
    tryall              : Try starting from each node(if time limit allows)\n\
"

static const char *nnOptionsStrings[] = {
    "random",
    "tryall"
};

#define ARGP_EM_DOC "\
Specify the way Extra Mileage heuristic is initialized every time it's called (DEFAULT=random)\n\
    random              : Use two random nodes\n\
    extremes            : Use a node with maximum x coordinate and one with minimum x coordinate\n\
    farthest            : Use farthest nodes\n\
"
//hull                : Compute the hull of the set using quickhull algorthim and use it as initialization.

static const char *emOptionsStrings[] = {
    "random",
    "extremes",
    "farthest"
};
// hull

#define ARGP_METAHEURISTICS_INIT_MODE_DOC "\
Specify the heuristic that vns shall use at the start when finding the base solution (DEFAULT=nn)\n\
    nn                  : Use Nearest Neighbor\n\
    em                  : Use Extra Mileage\n\
"
#define ARGP_METAHEURISTICS_INIT_MODES_COUNT HEURISTICS_MODES_COUNT

#define ARGP_MATHEUR_INIT_MODE_DOC "\
Specify which Heuristic/Metaheuristic to use as initialization for matheuristics\n\
    nn                  : Use Nearest Neighbor\n\
    em                  : Use Extra Milage\n\
    tabu                : Use Tabu Search\n\
    vns                 : Use Variable Neighborhood Search\n\
    genetic             : Use Genetic Heuristic\n\
    annealing           : Use Simulated Annealing\n\
"
#define ARGP_MATHEUR_INIT_MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT)

#define ARGP_WARM_START_MODE_DOC "\
Specify the type of heuristic to use when computing a warm-start solution\n\
    nn                  : Use Nearest Neighbor\n\
    em                  : Use Extra Milage\n\
    tabu                : Use Tabu Search\n\
    vns                 : Use Variable Neighborhood Search\n\
    genetic             : Use Genetic Heuristic\n\
    annealing           : Use Simulated Annealing\n\
"
#define ARGP_WARM_START_MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT)

#define ARGP_HARDFIX_POLICY_DOC "\
Specify the type of policy Hard Fixing should use when fixing edges\n\
    random              : Fix random edges\n\
    small               : Fix edges with smallest cost in solution\n\
    prob                : Fix edges with probability inversely proportional to their cost\n\
    mixed               : Use a mix of the policies above\n\
"

static const char *hardFixPolicyStrings[] = {
    "random",
    "small",
    "prob",
    "mixed"
};

#define LOG_LEVEL_DOC "\
Specify the log level (or verbosity level) of for the run. (DEFAULT=log\n\
    error   : Show only error messages\n\
    critical: Show critical messages and all above\n\
    warning : Show warning and all above\n\
    notice  : Show notice messages and all above\n\
    log     : Show log messages and all above\n\
    debug   : Show debug messages and all above\n\
    all     : Show all messages\n\
"

static const char *logLevelStrings[] = {
    "error",
    "critical",
    "warning",
    "notice",
    "log",
    "debug",
    "all"
};


enum argpKeys{
    ARGP_FILE='f',
    ARGP_MODE='m',

    ARGP_GRASP='g',
    ARGP_2OPT='2',
    ARGP_TLIM='t',

    ARGP_NN_MODE=257,
    ARGP_EM_MODE,

    ARGP_METAHEURISTICS_INIT_MODE,
    ARGP_MATHEURISTICS_INIT_MODE,
    ARGP_WARM_START_MODE,

    ARGP_HARDFIX_POLICY,

    ARGP_SEED,
    ARGP_NTHREADS='j',
    ARGP_ROUND='r',
    ARGP_PLOT='p',
    ARGP_SAVE='s',
    ARGP_LOG_LEVEL='l'
};

error_t argpParser(int key, char *arg, struct argp_state *state);

static int parseMode(char *arg, Instance *inst);

static int parseGrasp(char *arg, Instance *inst);
static int parseTlim(char *arg, Instance *inst);

static int parseNNOption(char *arg, Instance *inst);
static int parseEMOption(char *arg, Instance *inst);

static int parseMetaheurInitMode(char *arg, Instance *inst);
static int parseMatheurInitMode(char *arg, Instance *inst);
static int parseWarmStartMode(char *arg, Instance *inst);

static int parseHardFixPolicy(char *arg, Instance *inst);

static int parseSeed(char *arg, Instance *inst);
static int parseNThreads(char *arg, Instance *inst);
static int parseLogLevel(char *arg, Instance *inst);

static int checkEssentials(Instance *inst);


void argParse(Instance * inst, int argc, char *argv[])
{
    static struct argp_option argpOptions[] = {
        { .name="file", .key=ARGP_FILE, .arg="FILENAME", .flags=0, .doc="Location of the .tsp file containing the instance to use\n", .group=1 },
        { .name="mode", .key=ARGP_MODE, .arg="MODE", .flags=0, .doc=ARGP_MODE_DOC, .group=1 },

        { .name="grasp", .key=ARGP_GRASP, .arg="STRING", .flags=0, .doc=ARGP_GRASP_DOC, .group=2 },
        { .name="2opt", .key=ARGP_2OPT, .arg=NULL, .flags=0, .doc="Specify to use 2-opt at the end of the selected heuristic\n", .group=2 },
        { .name="tlim", .key=ARGP_TLIM, .arg="UINT", .flags=0, .doc="Specify time limit for the execution\n", .group=2 },

        { .name="nn-opt", .key=ARGP_NN_MODE, .arg="STRING", .flags=0, .doc=ARGP_NN_DOC, .group=3 },
        { .name="em-opt", .key=ARGP_EM_MODE, .arg="STRING", .flags=0, .doc=ARGP_EM_DOC, .group=3 },

        { .name="metaheur-init-mode", .key=ARGP_METAHEURISTICS_INIT_MODE, .arg="STRING", .flags=0, .doc=ARGP_METAHEURISTICS_INIT_MODE_DOC, .group=3 },
        { .name="matheur-init-mode", .key=ARGP_MATHEURISTICS_INIT_MODE, .arg="STRING", .flags=0, .doc=ARGP_MATHEUR_INIT_MODE_DOC, .group=3 },
        { .name="warm-start-mode", .key=ARGP_WARM_START_MODE, .arg="STRING", .flags=0, .doc=ARGP_WARM_START_MODE_DOC, .group=3 },

        { .name="hardfix-policy", .key=ARGP_HARDFIX_POLICY, .arg="STRING", .flags=0, .doc=ARGP_HARDFIX_POLICY_DOC, .group=3 },

        { .name="seed", .key=ARGP_SEED, .arg="UINT", .flags=0, .doc="Random Seed [0,MAX_INT32] to use as random seed for the current run. If -1 seed will be random\n", .group=4 },
        { .name="threads", .key=ARGP_NTHREADS, .arg="UINT", .flags=0, .doc="Maximum number of threads to use. If not specified gets maximum automatically\n", .group=4 },
        { .name="roundcosts", .key=ARGP_ROUND, .arg=NULL, .flags=0, .doc="Specify this if yout want to use rounded version of edge cost\n", .group=4 },
        { .name="plot", .key=ARGP_PLOT, .arg=NULL, .flags=0, .doc="Specify this if yout want to plot final result\n", .group=4 },
        { .name="save", .key=ARGP_SAVE, .arg=NULL, .flags=0, .doc="Specify this if yout want to save final result in run/\n", .group=4 },
        { .name="loglvl", .key=ARGP_LOG_LEVEL, .arg="STRING", .flags=0, .doc=LOG_LEVEL_DOC, .group=4 },
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
        return parseMode(arg, inst);

    case ARGP_GRASP:
        return parseGrasp(arg, inst);

    case ARGP_2OPT:
        inst->params.use2OptFlag = 1;
        break;
    
    case ARGP_TLIM:
        return parseTlim(arg, inst);

    case ARGP_NN_MODE:
        return parseNNOption(arg, inst);
    
    case ARGP_EM_MODE:
        return parseEMOption(arg, inst);
    
    case ARGP_METAHEURISTICS_INIT_MODE:
        return parseMetaheurInitMode(arg, inst);
    
    case ARGP_MATHEURISTICS_INIT_MODE:
        return parseMatheurInitMode(arg, inst);

    case ARGP_WARM_START_MODE:
        return parseWarmStartMode(arg, inst);

    case ARGP_HARDFIX_POLICY:
        return parseHardFixPolicy(arg, inst);

    case ARGP_SEED:
        return parseSeed(arg, inst);

    case ARGP_NTHREADS:
        return parseNThreads(arg, inst);

    case ARGP_ROUND:
        inst->params.roundWeightsFlag = 1;
        break;

    case ARGP_PLOT:
        inst->params.showPlotFlag = 1;
        break;
    
    case ARGP_SAVE:
        inst->params.saveFlag = 1;
        break;
    
    case ARGP_LOG_LEVEL:
        return parseLogLevel(arg, inst);
    
    case ARGP_KEY_END:
        // check if necessary flags have been provided
        return checkEssentials(inst);

    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static int parseMode(char *arg, Instance *inst)
{
    for (size_t i = 0; i < MODES_COUNT; i++)
    {
        if (strcmp(arg, modeStrings[i]) == 0)
        {
            inst->params.mode = i;
            return 0;
        }
    }

    LOG(LOG_LVL_ERROR, "Mode specified is not valid");
    return ARGP_ERR_UNKNOWN;
}

static int parseGrasp(char *arg, Instance *inst)
{
    if (strncmp(arg, graspStrings[GRASP_ALMOSTBEST], strlen(graspStrings[GRASP_ALMOSTBEST])) == 0)
    {
        inst->params.graspType = GRASP_ALMOSTBEST;

        // check if GRASP_CHANCE is also specified
        size_t graspChancePosStart = strlen(graspStrings[GRASP_ALMOSTBEST]);
        if (arg[graspChancePosStart] == '(')
        {
            char *endPtr;
            double cvt = strtod(&arg[graspChancePosStart + 1], &endPtr);
            if (*endPtr != ')')
            {
                LOG(LOG_LVL_ERROR, "Couldn't find the ')' after GRASP_CHANCE. It must be placed directly after the number");
                return ARGP_ERR_UNKNOWN;
            }
            if ((cvt <= 0) || (cvt > 1))
            {
                LOG(LOG_LVL_ERROR, "GRASP_CHANCE specified is not a valid number. It must be a number in (0,1]");
                return ARGP_ERR_UNKNOWN;
            }
            if (endPtr[1] != 0)
                LOG(LOG_LVL_WARNING, "There are elements after the ')'. They will be ignored");
            inst->params.graspChance = cvt;
        }
    }
    else if (strncmp(arg, graspStrings[GRASP_RANDOM], strlen(graspStrings[GRASP_RANDOM])) == 0)
    {
        inst->params.graspType = GRASP_RANDOM;

        // check if GRASP_CHANCE is also specified
        size_t graspChancePosStart = strlen(graspStrings[GRASP_RANDOM]);
        if (arg[graspChancePosStart] == '(')
        {
            char *endPtr;
            double cvt = strtod(&arg[graspChancePosStart + 1], &endPtr);
            if (*endPtr != ')')
            {
                LOG(LOG_LVL_ERROR, "Couldn't find the ')' after GRASP_CHANCE. It must be placed directly after the number");
                return ARGP_ERR_UNKNOWN;
            }
            if ((cvt <= 0) || (cvt > 1))
            {
                LOG(LOG_LVL_ERROR, "GRASP_CHANCE specified is not a valid number. It must be a number in (0,1]");
                return ARGP_ERR_UNKNOWN;
            }
            if (endPtr[1] != 0)
                LOG(LOG_LVL_WARNING, "There are elements after the ')'. They will be ignored");
            inst->params.graspChance = cvt;
        }
    }

    return 0;
}

static int parseTlim(char *arg, Instance *inst)
{
    char *endPtr;
    long cvt = strtol(arg, &endPtr, 10);
    if (cvt <= 0)
    {
        LOG(LOG_LVL_ERROR, "The number specified as time limit must be in seconds and be > 0");
        return ARGP_ERR_UNKNOWN;
    }
    if (endPtr != &arg[strlen(arg)])
        LOG(LOG_LVL_WARNING, "There are extra character after the time limit value");

    inst->params.tlim = (double)cvt;
    return 0;
}

static int parseNNOption(char *arg, Instance *inst)
{
    if (strcmp(arg, nnOptionsStrings[NN_FIRST_RANDOM]) == 0)
        inst->params.nnFirstNodeOption = NN_FIRST_RANDOM;
    else if (strcmp(arg, nnOptionsStrings[NN_FIRST_TRYALL]) == 0)
        inst->params.nnFirstNodeOption = NN_FIRST_TRYALL;
    else
    {
        LOG(LOG_LVL_ERROR, "nn-opt: argument not recognized");
        return ARGP_ERR_UNKNOWN;
    }
    
    return 0;
}

static int parseEMOption(char *arg, Instance *inst)
{
    if (strcmp(arg, emOptionsStrings[EM_INIT_RANDOM]) == 0)
        inst->params.emInitOption = EM_INIT_RANDOM;
    else if (strcmp(arg, emOptionsStrings[EM_INIT_EXTREMES]) == 0)
        inst->params.emInitOption = EM_INIT_EXTREMES;
    else if (strcmp(arg, emOptionsStrings[EM_INIT_FARTHEST_POINTS]) == 0)
        inst->params.emInitOption = EM_INIT_FARTHEST_POINTS;
    else
    {
        LOG(LOG_LVL_ERROR, "em-opt: argument not recognized");
        return ARGP_ERR_UNKNOWN;
    }
    
    return 0;
}

static int parseMetaheurInitMode(char *arg, Instance *inst)
{
    for (size_t i = 0; i < ARGP_METAHEURISTICS_INIT_MODES_COUNT; i++)
    {
        if (strcmp(arg, modeStrings[i]) == 0)
        {
            inst->params.metaheurInitMode = i;
            return 0;
        }
    }

    LOG(LOG_LVL_ERROR, "metaheur-init-mode: argument not recognized");
    return ARGP_ERR_UNKNOWN;
}

static int parseMatheurInitMode(char *arg, Instance *inst)
{
    for (size_t i = 0; i < ARGP_MATHEUR_INIT_MODES_COUNT; i++)
    {
        if (strcmp(arg, modeStrings[i]) == 0)
        {
            inst->params.matheurInitMode = i;
            return 0;
        }
    }

    LOG(LOG_LVL_ERROR, "math-init-mode: argument not recognized");
    return ARGP_ERR_UNKNOWN;
}


static int parseWarmStartMode(char *arg, Instance *inst)
{
    for (size_t i = 0; i < ARGP_WARM_START_MODES_COUNT; i++)
    {
        if (strcmp(arg, modeStrings[i]) == 0)
        {
            inst->params.warmStartMode = i;
            return 0;
        }
    }

    LOG(LOG_LVL_ERROR, "warm-start-mode: argument not recognized");
    return ARGP_ERR_UNKNOWN;
}

static int parseHardFixPolicy(char *arg, Instance *inst)
{
    if (strcmp(arg, hardFixPolicyStrings[HARDFIX_POLICY_RANDOM]) == 0)
        inst->params.hardFixPolicy = HARDFIX_POLICY_RANDOM;
    else if (strcmp(arg, hardFixPolicyStrings[HARDFIX_POLICY_SMALLEST]) == 0)
        inst->params.hardFixPolicy = HARDFIX_POLICY_SMALLEST;
    else if (strcmp(arg, hardFixPolicyStrings[HARDFIX_POLICY_PROBABILITY]) == 0)
        inst->params.hardFixPolicy = HARDFIX_POLICY_PROBABILITY;
    else if (strcmp(arg, hardFixPolicyStrings[HARDFIX_POLICY_MIXED]) == 0)
        inst->params.hardFixPolicy = HARDFIX_POLICY_MIXED;
    else
    {
        LOG(LOG_LVL_ERROR, "hardfix-policy: argument not recognized");
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static int parseSeed(char *arg, Instance *inst)
{
    char *endPtr;
    long cvt = strtol(arg, &endPtr, 10);
    if (cvt <= 0)
    {
        LOG(LOG_LVL_ERROR, "The number specified as seed must be an interger value > 0");
        return ARGP_ERR_UNKNOWN;
    }
    if (endPtr != &arg[strlen(arg)])
        LOG(LOG_LVL_WARNING, "There are extra character after the seed value");

    inst->params.randomSeed = (int)cvt;
    return 0;
}

static int parseNThreads(char *arg, Instance *inst)
{
    char *endPtr;
    long cvt = strtol(arg, &endPtr, 10);
    if (cvt <= 0 || cvt > 32)
    {
        LOG(LOG_LVL_ERROR, "The number of threads must be an interger value > 0 and < 32");
        return ARGP_ERR_UNKNOWN;
    }
    if (endPtr != &arg[strlen(arg)])
        LOG(LOG_LVL_WARNING, "There are extra character after the thread value");

    inst->params.nThreads = (int)cvt;
    return 0;
}

static int parseLogLevel(char *arg, Instance *inst)
{
    if (strcmp(arg, logLevelStrings[LOG_LVL_ERROR]) == 0)
        inst->params.logLevel = LOG_LVL_ERROR;
    else if (strcmp(arg, logLevelStrings[LOG_LVL_CRITICAL]) == 0)
        inst->params.logLevel = LOG_LVL_CRITICAL;
    else if (strcmp(arg, logLevelStrings[LOG_LVL_WARNING]) == 0)
        inst->params.logLevel = LOG_LVL_WARNING;
    else if (strcmp(arg, logLevelStrings[LOG_LVL_NOTICE]) == 0)
        inst->params.logLevel = LOG_LVL_NOTICE;
    else if (strcmp(arg, logLevelStrings[LOG_LVL_LOG]) == 0)
        inst->params.logLevel = LOG_LVL_LOG;
    else if (strcmp(arg, logLevelStrings[LOG_LVL_DEBUG]) == 0)
        inst->params.logLevel = LOG_LVL_DEBUG;
    else if (strcmp(arg, logLevelStrings[LOG_LVL_EVERYTHING]) == 0)
        inst->params.logLevel = LOG_LVL_EVERYTHING;
    else
    {
        LOG(LOG_LVL_ERROR, "loglvl: argument not recognized");
        return ARGP_ERR_UNKNOWN;
    }

    setLogLevel(inst->params.logLevel);

    return 0;
}

static int checkEssentials(Instance *inst)
{
    if (inst->params.mode == MODE_NONE)
    {
        LOG(LOG_LVL_ERROR, "Missing --mode (-m) option in arguments. A mode must be specified. Check --help or --usage to check what modes are avilable");
        return ARGP_ERR_UNKNOWN;
    }
    if (inst->params.inputFile[0] == 0)
    {
        LOG(LOG_LVL_ERROR, "Missing --file (-f) option in arguments. An input file must be specified");
        return ARGP_ERR_UNKNOWN;
    }
    if (inst->params.tlim == -1)
    {
        if (inst->params.graspType != GRASP_NONE)
        {
            LOG(LOG_LVL_ERROR, "If grasp is used a time limit must be provided with option --tlim (-t)");
            return ARGP_ERR_UNKNOWN;
        }
        if (inst->params.mode >= MODE_VNS)
        {
            LOG(LOG_LVL_ERROR, "When using an heuristic that is intrinsically using grasp(like %s), a time limit must also be specified", modeStrings[inst->params.mode]);
            return ARGP_ERR_UNKNOWN;
        }
    }
    return 0;
}


void printInfo(Instance *inst)
{
    const char * blank = "   ";

    Parameters *p = &inst->params;

    printf("SETTINGS:\n");

    // input file
    printf("%sInput File/Problem: \"%s\"\n", blank, p->inputFile);
    // mode
    printf("%sCurrent running mode is %s\n", blank, modeStrings[p->mode]);

    // grasp
    if (p->graspType == GRASP_NONE)
        printf("%sGrasp is off\n", blank);
    else
        printf("%sGrasp is on, grasp mode id is %s with chance %lf\n", blank, graspStrings[p->graspType], p->graspChance);
    // 2Opt
    if (p->use2OptFlag || p->mode == MODE_VNS)
        printf("%sUsing 2Opt\n", blank);
    // time limit
    if (p->tlim != -1)
        printf("%sTime limit of %lf seconds\n", blank, p->tlim);
    else
        printf("%sTime limit is not set\n", blank);
    
    // nn options
    if (p->mode == MODE_NN || (p->mode >= HEURISTICS_MODES_COUNT  && (p->metaheurInitMode == MODE_NN && p->matheurInitMode == MODE_EM)) || 
            (p->mode >= (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT) && p->matheurInitMode == MODE_NN))
        printf("%sNearest Neighbor starting node set to: %s\n", blank, nnOptionsStrings[p->nnFirstNodeOption]);
    // em options
    if (p->mode == MODE_EM || (p->mode >= HEURISTICS_MODES_COUNT  && (p->metaheurInitMode == MODE_EM && p->matheurInitMode == MODE_NN)) || 
            (p->mode >= (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT) && p->matheurInitMode == MODE_EM))
        printf("%sExtra Mileage initialization set to: %s\n", blank, emOptionsStrings[p->emInitOption]);
    // metaheuristics modes
    if (p->mode >= HEURISTICS_MODES_COUNT && (p->mode < HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT || p->matheurInitMode >= HEURISTICS_MODES_COUNT))
        printf("%sMetaheuristics initialization set to: %s\n", blank, modeStrings[p->metaheurInitMode]);
    // matheuristics modes
    if (p->mode >= HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT)
        printf("%sMatheuristics initialization set to: %s\n", blank, modeStrings[p->matheurInitMode]);

    // seed
    if (p->randomSeed != -1)
        printf("%sRandom Seed = %d\n", blank, p->randomSeed);
    // threads
    printf("%sThreads used = %d\n", blank, p->nThreads);
    // roundcosts
    if (p->roundWeightsFlag)
        printf("%sEdge Cost is rounded according to tsplib documentation file\n", blank);
    // save
    if (p->saveFlag)
        printf("%sFinal solution of this run will be saved in a .tour file inside OperationsResearch2/runs\n", blank);
    // log level
    printf("%sLog level = %s", blank, logLevelStrings[p->logLevel]);

    printf("\n");
}
