#include "Tsp.h"

#include <argp.h>
#include <unistd.h>
#include <string.h>

#define SUBOPT_BLANKSPACE "  "

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


#define DOC_NN SUBOPT_BLANKSPACE SUBOPT_NN "          : Use Nearest Neighbor\n"
#define DOC_EM SUBOPT_BLANKSPACE SUBOPT_EM "          : Use Extra Mileage\n"
#define DOC_TABU SUBOPT_BLANKSPACE SUBOPT_TABU "        : Use Tabu Search\n"
#define DOC_VNS SUBOPT_BLANKSPACE SUBOPT_VNS "         : Use Variable Neighborhood Search\n"
#define DOC_ANNEALING SUBOPT_BLANKSPACE SUBOPT_ANNEALING "   : Use Simulated Annealing\n"
#define DOC_GENETIC SUBOPT_BLANKSPACE SUBOPT_GENETIC "     : Use Genetic Heuristic\n"
#define DOC_BENDERS SUBOPT_BLANKSPACE SUBOPT_BENDERS "     : Use Benders method\n"
#define DOC_BRANCHCUT SUBOPT_BLANKSPACE SUBOPT_BRANCHCUT "  : Use cplex generic callback to perform Branch & Cut\n"
#define DOC_HARDFIX SUBOPT_BLANKSPACE SUBOPT_HARDFIX "     : Use Hard Fixing Matheuristic\n"
#define DOC_LOCALBRANCHING SUBOPT_BLANKSPACE SUBOPT_LOCALBRANCHING " : Use Local Branching Matheuristic\n"

#define ARGP_MODE_DOC "\
Specify the type of solver to use: \n"\
DOC_NN DOC_EM DOC_TABU DOC_VNS DOC_ANNEALING DOC_GENETIC DOC_BENDERS DOC_BRANCHCUT DOC_HARDFIX DOC_LOCALBRANCHING

#define HEURISTICS_MODES_COUNT 2
#define METAHEUR_MODES_COUNT 4
#define EXACT_SOLVERS_COUNT 2
#define MATHEURISTICS_COUNT 2

#define MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT + EXACT_SOLVERS_COUNT + MATHEURISTICS_COUNT)
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
Specify to Grasp mode (DEFAULT=random(0.1))\n" \
SUBOPT_BLANKSPACE SUBOPT_GRASP_ALMOSTBEST " : The Use of grasp will be limited to selecting another good choice with default probability value\n" \
SUBOPT_BLANKSPACE SUBOPT_GRASP_ALMOSTBEST "[<GRASP_CHANCE>] : Same as \"almostbest\" but the selection of the secondary choice happens with probability specified with <GRASP_CHANCE>\n" \
SUBOPT_BLANKSPACE SUBOPT_GRASP_RANDOM "     : At every iteration have a completely random choice with default probability\n" \
SUBOPT_BLANKSPACE SUBOPT_GRASP_RANDOM "[<GRASP_CHANCE>]     : Same as \"random\" but the probability of the random choice is specified with <GRASP_CHANCE>\n"

static const char *graspStrings[] = { SUBOPT_GRASP_ALMOSTBEST, SUBOPT_GRASP_RANDOM };

#define METAHEURISTICS_INIT_MODE_DOC "\
Specify the heuristic that vns shall use at the start when finding the base solution (DEFAULT=nn)\n" \
DOC_NN DOC_EM
#define METAHEURISTICS_INIT_MODES_COUNT HEURISTICS_MODES_COUNT

#define MATHEUR_INIT_MODE_DOC "\
Specify which Heuristic/Metaheuristic to use as initialization for matheuristics\n" \
DOC_NN DOC_EM DOC_TABU DOC_VNS DOC_ANNEALING DOC_GENETIC
#define MATHEUR_INIT_MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT)

#define WARM_START_MODE_DOC "\
Specify the type of heuristic to use when computing a warm-start solution\n" \
DOC_NN DOC_EM DOC_TABU DOC_VNS DOC_ANNEALING DOC_GENETIC
#define WARM_START_MODES_COUNT (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT)

#define LOG_LEVEL_DOC "\
Specify the log level (or verbosity level) of for the run. (DEFAULT=log\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_ERROR "   : Show only error messages\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_CRITICAL ": Show critical messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_WARNING " : Show warning and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_NOTICE "  : Show notice messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_LOG "     : Show log messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_DEBUG "   : Show debug messages and all above\n" \
SUBOPT_BLANKSPACE SUBOPT_LOG_ALL "     : Show all messages\n"

static const char *logLevelStrings[] = { SUBOPT_LOG_ERROR, SUBOPT_LOG_CRITICAL, SUBOPT_LOG_WARNING, SUBOPT_LOG_NOTICE, SUBOPT_LOG_LOG, SUBOPT_LOG_DEBUG, SUBOPT_LOG_ALL };
static const int loglvlsCount = sizeof(logLevelStrings)/sizeof(*logLevelStrings);

enum argpKeys{
    ARGP_FILE='f',
    ARGP_MODE='m',

    ARGP_GRASP='g',
    ARGP_2OPT='2',
    ARGP_TLIM='t',

    ARGP_METAHEURISTICS_INIT_MODE=257,
    ARGP_MATHEURISTICS_INIT_MODE,
    ARGP_WARM_START_MODE,

    ARGP_NN_TRYALL,
    ARGP_EM_FARTHEST,
    ARGP_HARDFIX_SMALLEST,

    ARGP_SEED,
    ARGP_NTHREADS='j',
    ARGP_ROUND='r',
    ARGP_PLOT='p',
    ARGP_SAVE='s',
    ARGP_LOG_LEVEL='l'
};

error_t argpParser(int key, char *arg, struct argp_state *state);

static int parseEnumOption(char *arg, int *savePtr, const char **optionsSet, const int from, const int to, const char *optionName);

static int parseGrasp(char *arg, Instance *inst);
static int parseTlim(char *arg, Instance *inst);

static int parseSeed(char *arg, Instance *inst);
static int parseNThreads(char *arg, Instance *inst);

static int checkEssentials(Instance *inst);


void argParse(Instance * inst, int argc, char *argv[])
{
    static struct argp_option argpOptions[] = {
        { .name="file", .key=ARGP_FILE, .arg="FILENAME", .flags=0, .doc="Location of the .tsp file containing the instance to use\n", .group=1 },
        { .name="mode", .key=ARGP_MODE, .arg="MODE", .flags=0, .doc=ARGP_MODE_DOC, .group=1 },

        { .name="grasp", .key=ARGP_GRASP, .arg="STRING", .flags=0, .doc=GRASP_DOC, .group=2 },
        { .name="2opt", .key=ARGP_2OPT, .arg=NULL, .flags=0, .doc="Specify to use 2-opt at the end of the selected heuristic\n", .group=2 },
        { .name="tlim", .key=ARGP_TLIM, .arg="UINT", .flags=0, .doc="Specify time limit for the execution\n", .group=2 },

        { .name="metaheur-init-mode", .key=ARGP_METAHEURISTICS_INIT_MODE, .arg="STRING", .flags=0, .doc=METAHEURISTICS_INIT_MODE_DOC, .group=3 },
        { .name="matheur-init-mode", .key=ARGP_MATHEURISTICS_INIT_MODE, .arg="STRING", .flags=0, .doc=MATHEUR_INIT_MODE_DOC, .group=3 },
        { .name="warm-start-mode", .key=ARGP_WARM_START_MODE, .arg="STRING", .flags=0, .doc=WARM_START_MODE_DOC, .group=3 },

        { .name="nn-tryall", .key=ARGP_NN_TRYALL, .arg=NULL, .flags=0, .doc="Specify to make Nearest Neighbor start from each node instead of chosing a random one each time\n", .group=4 },
        { .name="em-farthest", .key=ARGP_EM_FARTHEST, .arg=NULL, .flags=0, .doc="Specify to make Extra Mileage initialization the farthest nodes each time instead of a random one each time\n", .group=4 },
        { .name="hardfix-smallest", .key=ARGP_HARDFIX_SMALLEST, .arg=NULL, .flags=0, .doc="Specify to make Hard Fixing fix only edges with smallest cost instead of fixing random edges\n", .group=4 },

        { .name="seed", .key=ARGP_SEED, .arg="UINT", .flags=0, .doc="Random Seed [0,MAX_INT32] to use as random seed for the current run. If -1 seed will be random\n", .group=5 },
        { .name="threads", .key=ARGP_NTHREADS, .arg="UINT", .flags=0, .doc="Maximum number of threads to use. If not specified gets maximum automatically\n", .group=5 },
        { .name="roundcosts", .key=ARGP_ROUND, .arg=NULL, .flags=0, .doc="Specify this if yout want to use rounded version of edge cost\n", .group=5 },
        { .name="plot", .key=ARGP_PLOT, .arg=NULL, .flags=0, .doc="Specify this if yout want to plot final result\n", .group=5 },
        { .name="save", .key=ARGP_SAVE, .arg=NULL, .flags=0, .doc="Specify this if yout want to save final result in run/\n", .group=5 },
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
    int retVal = 0;

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
        return parseEnumOption(arg, &inst->params.mode, modeStrings, 0, MODES_COUNT, "mode");

    case ARGP_GRASP:
        return parseGrasp(arg, inst);

    case ARGP_2OPT:
        inst->params.use2Opt = true;
        break;
    
    case ARGP_TLIM:
        return parseTlim(arg, inst);
    
    case ARGP_METAHEURISTICS_INIT_MODE:
        return parseEnumOption(arg, (int*)&inst->params.metaheurInitMode, modeStrings, 0, HEURISTICS_MODES_COUNT, "metaheur-init-mode");
    
    case ARGP_MATHEURISTICS_INIT_MODE:
        return parseEnumOption(arg, (int*)&inst->params.matheurInitMode, modeStrings, 0, HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT, "matheur-init-mode");

    case ARGP_WARM_START_MODE:
        return parseEnumOption(arg, (int*)&inst->params.warmStartMode, modeStrings, 0, HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT, "warm-start-mode");

    case ARGP_NN_TRYALL:
        inst->params.nnFirstNodeOption = NN_FIRST_TRYALL;
        break;
    
    case ARGP_EM_FARTHEST:
        inst->params.emInitOption = EM_INIT_FARTHEST_POINTS;
        break;

    case ARGP_HARDFIX_SMALLEST:
        inst->params.hardFixPolicy = HARDFIX_POLICY_SMALLEST;
        break;

    case ARGP_SEED:
        return parseSeed(arg, inst);

    case ARGP_NTHREADS:
        return parseNThreads(arg, inst);

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
        retVal = parseEnumOption(arg, (int*)&inst->params.logLevel, logLevelStrings, 0, loglvlsCount, "loglvl");
        setLogLevel(inst->params.logLevel);
        return retVal;
    
    case ARGP_KEY_END:
        // check if necessary flags have been provided
        return checkEssentials(inst);

    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static int parseEnumOption(char *arg, int *savePtr, const char **optionsSet, const int from, const int to, const char *optionName)
{
    for (int i = from; i < to; i++)
    {
        if (strcmp(arg, optionsSet[i]) == 0)
        {
            *savePtr = i;
            return 0;
        }
    }
    
    LOG(LOG_LVL_ERROR, "%s: argument not valid", optionName);
    return ARGP_ERR_UNKNOWN;
}

static int parseGrasp(char *arg, Instance *inst)
{
    if (strncmp(arg, graspStrings[GRASP_ALMOSTBEST], strlen(graspStrings[GRASP_ALMOSTBEST])) == 0)
    {
        inst->params.graspType = GRASP_ALMOSTBEST;

        // check if GRASP_CHANCE is also specified
        int graspChancePosStart = (int)strlen(graspStrings[GRASP_ALMOSTBEST]);
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
        int graspChancePosStart = (int)strlen(graspStrings[GRASP_RANDOM]);
        if (arg[graspChancePosStart] == '[')
        {
            char *endPtr;
            double cvt = strtod(&arg[graspChancePosStart + 1], &endPtr);
            if (*endPtr != ']')
            {
                LOG(LOG_LVL_ERROR, "Couldn't find the ']' after GRASP_CHANCE. It must be placed directly after the number");
                return ARGP_ERR_UNKNOWN;
            }
            if ((cvt <= 0) || (cvt > 1))
            {
                LOG(LOG_LVL_ERROR, "GRASP_CHANCE specified is not a valid number. It must be a number in (0,1]");
                return ARGP_ERR_UNKNOWN;
            }
            if (endPtr[1] != 0)
                LOG(LOG_LVL_WARNING, "There are elements after the ']'. They will be ignored");
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
        LOG(LOG_LVL_ERROR, "The number of threads must be an interger value > 0 and < %d", MAX_THREADS);
        return ARGP_ERR_UNKNOWN;
    }
    if (endPtr != &arg[strlen(arg)])
        LOG(LOG_LVL_WARNING, "There are extra character after the thread value");

    inst->params.nThreads = (int)cvt;
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


#define BLANK_SPACE "   "
void printInfo(Instance *inst)
{
    Parameters *p = &inst->params;

    printf("SETTINGS:\n");

    // input file
    printf(BLANK_SPACE "Input File/Problem: \"%s\"\n", p->inputFile);
    // mode
    printf(BLANK_SPACE "Current running mode is %s\n", modeStrings[p->mode]);

    // grasp
    if (p->graspType == GRASP_NONE)
        printf(BLANK_SPACE "Grasp is off\n");
    else
        printf(BLANK_SPACE "Grasp is on, grasp mode id is %s with chance %lf\n", graspStrings[p->graspType], p->graspChance);
    // 2Opt
    if (p->use2Opt || p->mode == MODE_VNS)
        printf(BLANK_SPACE "Using 2Opt\n");
    // time limit
    if (p->tlim != -1)
        printf(BLANK_SPACE "Time limit of %lf seconds\n", p->tlim);
    else
        printf(BLANK_SPACE "Time limit is not set\n");
    
    // metaheuristics modes
    if (p->mode >= HEURISTICS_MODES_COUNT && (p->mode < HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT || p->matheurInitMode >= HEURISTICS_MODES_COUNT))
        printf(BLANK_SPACE "Metaheuristics initialization set to: %s\n", modeStrings[p->metaheurInitMode]);
    // matheuristics modes
    if (p->mode >= HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT)
        printf(BLANK_SPACE "Matheuristics initialization set to: %s\n", modeStrings[p->matheurInitMode]);

    // nn options
    if (p->mode == MODE_NN || (p->mode >= HEURISTICS_MODES_COUNT  && (p->metaheurInitMode == MODE_NN && p->matheurInitMode == MODE_EM)) || 
            (p->mode >= (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT) && p->matheurInitMode == MODE_NN))
        printf(BLANK_SPACE "Nearest Neighbor starting node set to: %s\n", inst->params.nnFirstNodeOption == NN_FIRST_RANDOM ? "random" : "tryall");
    // em options
    if (p->mode == MODE_EM || (p->mode >= HEURISTICS_MODES_COUNT  && (p->metaheurInitMode == MODE_EM && p->matheurInitMode == MODE_NN)) || 
            (p->mode >= (HEURISTICS_MODES_COUNT + METAHEUR_MODES_COUNT) && p->matheurInitMode == MODE_EM))
        printf(BLANK_SPACE "Extra Mileage initialization set to: %s\n", inst->params.emInitOption == EM_INIT_RANDOM ? "random" : "farthest");
    // hardfix options
    if (p->mode == MODE_HARDFIX)
        printf(BLANK_SPACE "Hard Fixing policy set to: %s\n", inst->params.hardFixPolicy == HARDFIX_POLICY_RANDOM ? "random" : "smallest");

    // seed
    if (p->randomSeed != -1)
        printf(BLANK_SPACE "Random Seed = %d\n", p->randomSeed);
    // threads
    printf(BLANK_SPACE "Threads used = %d\n", p->nThreads);
    // roundcosts
    if (p->roundWeights)
        printf(BLANK_SPACE "Edge Cost is rounded according to tsplib documentation file\n");
    // save
    if (p->saveSolution)
        printf(BLANK_SPACE "Final solution of this run will be saved in a .tour file inside OperationsResearch2/runs\n");
    // log level
    printf(BLANK_SPACE "Log level = %s", logLevelStrings[p->logLevel]);

    printf("\n");
}
