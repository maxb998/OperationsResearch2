#include "ArgParser.h"
#include "TspUtilities.h"

#include <argp.h>
#include <unistd.h>
#include <string.h>


error_t argpParser(int key, char *arg, struct argp_state *state);

#define ARGP_MODE_DOC \
"Specify the type of solver to use: \n\
    nn  : Use Nearest Neighbor\n\
    em  : Use Extra Mileage\n\
    vns : Use Variable Neighborhood Search\
"
# define MODES_COUNT 3
const char *modeStrings[] = {
    "nn",
    "em",
    "vns"
};

#define ARGP_GRASP_DOC \
"Specify to use Grasp in the selected mode. Also set at least one Grasp setting between the followings\
    almostbest         : The Use of grasp will be limited(if possible) to selecting another good choice with default probability value \n\
    almostbest(<GRASP_CHANCE>) : Same as \"almostbest\" but the selection of the secondary choice happens with probability specified with <GRASP_CHANCE>\n\
    random             : At every iteration have a completely random choice with default probability\n\
    random(<GRASP_CHANCE>)     : Same as \"random\" but the probability of the random choice is specified with <GRASP_CHANCE>\n\
"

const char *graspStrings[] = {
    "almostbest",
    "random"
};


enum argpKeys{
    ARGP_FILE='f',
    ARGP_MODE='m',

    ARGP_GRASP='g',
    ARGP_2OPT='2',
    ARGP_TLIM='t',

    ARGP_SEED=256,
    ARGP_NTHREADS='j',
    ARGP_ROUND='r',
    ARGP_PLOT='p',
    ARGP_SAVE='s'
};


void argParse(Instance * inst, int argc, char *argv[])
{
    static struct argp_option argpOptions[] = {
        { .name="file", .key=ARGP_FILE, .arg="FILENAME", .flags=0, .doc="Location of the .tsp file containing the instance to use", .group=1 },
        { .name="mode", .key=ARGP_MODE, .arg="MODE", .flags=0, .doc=ARGP_MODE_DOC, .group=1 },

        { .name="grasp", .key=ARGP_GRASP, .arg="GRASP_TYPE", .flags=0, .doc=ARGP_GRASP_DOC, .group=2 },
        { .name="2opt", .key=ARGP_2OPT, .arg=NULL, .flags=0, .doc="Specify to use 2-opt at the end of the selected heuristic", .group=2 },
        { .name="tlim", .key=ARGP_TLIM, .arg="SECONDS", .flags=0, .doc="Specify time limit for the execution", .group=2 },

        { .name="seed", .key=ARGP_SEED, .arg="SEED", .flags=0, .doc="Random Seed [0,MAX_INT32] to use as random seed for the current run. If -1 seed will be random", .group=3 },
        { .name="threads", .key=ARGP_NTHREADS, .arg="N_THREADS", .flags=0, .doc="Maximum number of threads to use. If not specified gets maximum automatically", .group=3 },
        { .name="round", .key=ARGP_ROUND, .arg=NULL, .flags=0, .doc="Specify this if yout want to use rounded version of edge cost", .group=3 },
        { .name="plot", .key=ARGP_PLOT, .arg=NULL, .flags=0, .doc="Specify this if yout want to plot final result", .group=3 },
        { .name="save", .key=ARGP_SAVE, .arg=NULL, .flags=0, .doc="Specify this if yout want to save final result in run/", .group=3 },
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
        inst->params.mode = MODE_NONE;

        if (strcmp(arg, modeStrings[MODE_NN]) == 0)
            inst->params.mode = MODE_NN;
        else if (strcmp(arg, modeStrings[MODE_EM]) == 0)
            inst->params.mode = MODE_EM;
        else if (strcmp(arg, modeStrings[MODE_VNS]) == 0)
            inst->params.mode = MODE_VNS;
        
        // error check
        if (inst->params.mode == MODE_NONE)
        {
            LOG(LOG_LVL_ERROR, "Mode specified is not valid");
            return ARGP_ERR_UNKNOWN;
        }

    case ARGP_GRASP:
        if (strncmp(arg, graspStrings[GRASP_ALMOSTBEST], strlen(graspStrings[GRASP_ALMOSTBEST])) == 0)
        {
            inst->params.graspType = GRASP_ALMOSTBEST;

            // check if GRASP_CHANCE is also specified
            size_t graspChancePosStart = strlen(graspStrings[GRASP_ALMOSTBEST]);
            if (arg[graspChancePosStart] == '(')
            {
                char *endPtr;
                double cvt = strtod(&arg[graspChancePosStart+1], &endPtr);
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
                double cvt = strtod(&arg[graspChancePosStart+1], &endPtr);
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
        
        break;

    case ARGP_2OPT:
        inst->params.use2OptFlag = 1;
        break;
    
    case ARGP_TLIM:
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
        
        inst->params.tlim = (int) cvt;

        break;
    }

    case ARGP_SEED:
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
        
        inst->params.randomSeed = (int) cvt;
        break;
    }

    case ARGP_NTHREADS:
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
        
        inst->params.nThreads = (int) cvt;
        break;
    }

    case ARGP_ROUND:
        inst->params.roundWeightsFlag = 1;
        break;

    case ARGP_PLOT:
        inst->params.showPlotFlag = 1;
        break;
    
    case ARGP_SAVE:
        inst->params.saveFlag = 1;
        break;
    
    default:
        return ARGP_ERR_UNKNOWN;
    }

    return 0;
}