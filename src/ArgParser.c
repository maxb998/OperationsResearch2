#include "ArgParser.h"

#include <argp.h>

enum argpKeys{
    ARGP_FILE='f',
    ARGP_SEED='s',
    ARGP_NTHREADS='j',
    ARGP_ROUND='r',
    ARGP_PLOT='p'
};


void argParse(Instance * inst, int argc, char *argv[])
{
    static struct argp_option argpOptions[] = {
        { .name="file", .key=ARGP_FILE, .arg="FILENAME", .flags=0, .doc="Location of the .tsp file containing the instance to use", .group=0 },
        { .name="seed", .key=ARGP_SEED, .arg="SEED", .flags=0, .doc="Random Seed [0,MAX_INT32] to use as random seed for the current run. If -1 seed will be random", .group=1 },
        { .name="threads", .key=ARGP_NTHREADS, .arg="N_THREADS", .flags=0, .doc="Maximum number of threads to use. If not specified gets maximum automatically", .group=1 },
        { .name="round", .key=ARGP_ROUND, .arg=NULL, .flags=0, .doc="Specify this if yout want to use rounded version of edge cost", .group=1 },
        { .name="plot", .key=ARGP_PLOT, .arg=NULL, .flags=0, .doc="Specify this if yout want to plot final result", .group=1 },
        { .name="round", .key=ARGP_ROUND, .arg=NULL, .flags=0, .doc="Specify this if yout want to use rounded version of edge cost", .group=1 },
    };

    static struct argp argpData = {

    };
}