#ifndef UTILITIES
#define UTILITIES

#include "TspBase.h"

enum logLevel{
	LOG_LVL_ERROR, // 0
	LOG_LVL_CRITICAL, // 1
	LOG_LVL_WARNING, // 2
	LOG_LVL_NOTICE, // 3
	LOG_LVL_LOG, // 4
	LOG_LVL_DEBUG, // 5
	LOG_LVL_EVERYTHING // 6
};

// Swap elem1 and elem2. Can be done with any type of variable, however a temporary variable "tmp" of the same type of elem1 and elem2 MUST be provided.
#define swapElems(elem1,elem2,tmp) tmp = elem1; elem1 = elem2; elem2 = tmp

// Returns initialized/empty instance
Instance newInstance ();

// Returns initialized solution struct for the specified Instance pointed by inst
Solution newSolution (Instance *inst);

// Frees the allocated memory in the Instance pointed by d
void destroyInstance (Instance *inst);

// Frees the allocated memory in the Solution pointed by s
void destroySolution (Solution *sol);

Solution cloneSolution(Solution *sol);

/*Function used to log informations with level.
* lvl -> logLevel: Desired level of logging priority. If this value is greater than global value than LOG does not print anything
* line, ... -> string, args: Arguments that are directly fed to printf
*/
int LOG (enum logLevel lvl, char * line, ...);

// Launch fatal error, free memory and exit with code 1
void throwError (Instance *inst, Solution *sol, char * line, ...);

// Parse commandline arguments stored in argv and save relevant information to d->params
void parseArgs (Instance *inst, int argc, char *argv[]);

/* Checks on sol.indexPath that:
	1. First and last node in the path are the same (closed circuit)
	2. Nodes are not getting covered more than once
	3. All nodes are covered along the path*/
void basicSolutionCheck(Solution *sol);

/* Checks that the solution is feasible.
	Checks sol.X, sol.Y, sol.indexPath and also that the cost is correct*/
void fullSolutionCheck(Solution *sol);

// Check the correctness of the cost of the solution stored in Instance inst.
void costCheck(Solution *sol);

/*Plot solution using gnuplot. Does NOT check for errors on input
* d	-> Instance to plot
* plotPixelSize	-> string: Plot window size in pixel specified with format: "<WIDTH>,<HEIGHT>"
* pointColor -> string: Color of the circle representing the point, eg "black" or "red"
* tourPointColor -> string: Color of the 'X' on top of the point circle of color pointColor. Format and types is the same as for pointColor
* pointSize -> int: Size of the points
*/
void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize);

float computeSolutionCostVectorizedFloat(Solution *sol);

double computeSolutionCostVectorizedDouble(Solution *sol);

double computeSolutionCost(Solution *sol);

#endif // UTILITIES