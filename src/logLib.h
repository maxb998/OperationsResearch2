#ifndef LOGGER
#define LOGGER

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
 
enum logLevel{
	LOG_LVL_NONE, // 0
	LOG_LVL_CRITICAL, // 1
	LOG_LVL_WARNING, // 2
	LOG_LVL_NOTICE, // 3
	LOG_LVL_LOG, // 4
	LOG_LVL_DEBUG, // 5
	LOG_LVL_EVERYTHING // 6
};

const char * log_level_strings [] = {
	"NONE", // 0
	"CRIT", // 1
	"WARN", // 2
	"NOTI", // 3
	"LOG ", // 4
	"DEBG" // 5
};



int LOG (enum logLevel lvl, char * line, ...);

#endif //LOGGER