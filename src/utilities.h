#ifndef LOGGER
#define LOGGER

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <time.h>
#include <string.h>
#include <getopt.h>
 
#include "tsp.h"

enum logLevel{
	LOG_LVL_NONE, // 0
	LOG_LVL_CRITICAL, // 1
	LOG_LVL_WARNING, // 2
	LOG_LVL_NOTICE, // 3
	LOG_LVL_LOG, // 4
	LOG_LVL_DEBUG, // 5
	LOG_LVL_EVERYTHING // 6
};

int LOG (enum logLevel lvl, char * line, ...);

void parseArgs (Instance *d, int argc, char *argv[]);

void readFile (Instance *d);

#endif //LOGGER