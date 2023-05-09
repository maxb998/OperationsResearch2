#ifndef TSP_FUNCTIONS_H
#define TSP_FUNCTIONS_H

#include "TspBase.h"
#include "EdgeCostFunctions.h"

//###################################################################################################################################
// ARG_PARSER 
//###################################################################################################################################

void argParse(Instance * inst, int argc, char *argv[]);

void printInfo(Instance *inst);

//###################################################################################################################################
// TSP_IOUTILS
//###################################################################################################################################

// Read file with .tsp extension according to tsplib specifications, complete with file sintax error checking.
// Returns time elapsed while reading file
double readFile (Instance *inst);

void saveSolution(Solution *sol, int argc, char *argv[]);

/*Plot solution using gnuplot. Does NOT check for errors on input
 * d	-> Instance to plot
 * plotPixelSize	-> string: Plot window size in pixel specified with format: "<WIDTH>,<HEIGHT>"
 * pointColor -> string: Color of the circle representing the point, eg "black" or "red"
 * tourPointColor -> string: Color of the 'X' on top of the point circle of color pointColor. Format and types is the same as for pointColor
 * pointSize -> int: Size of the points
 * printIndex -> int: set to 1 to print index of each point as label on plot
*/
void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize, const int printIndex);

//###################################################################################################################################
// UTILITIES
//###################################################################################################################################

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

// Creates a copy of the solution passed. Coordinates and other parts that rely on malloc are completely copied element by element
Solution cloneSolution(Solution *sol);

// Set the log level for the LOG function. Must be called to set a log level that is different than LOG_LVL_LOG
void setLogLevel(enum logLevel lvl);

/*Function used to log informations with level.
* lvl -> logLevel: Desired level of logging priority. If this value is greater than global value than LOG does not print anything
* line, ... -> string, args: Arguments that are directly fed to printf
*/
int LOG (enum logLevel lvl, char * line, ...);

// Launch fatal error, destroys inst and sol (if not NULL) and exit with code 1
void throwError (Instance *inst, Solution *sol, char * line, ...);

void checkSolution(Solution *sol);

double computeSolutionCostVectorized(Solution *sol);

double computeSolutionCost(Solution *sol);


//###################################################################################################################################
// COST_MATRIX
//###################################################################################################################################

void printCostMatrix(Instance *inst);

double computeCostMatrix(Instance *inst);


//###################################################################################################################################
// NEAREST_NEIGHBOR
//###################################################################################################################################

Solution NearestNeighbor(Instance *inst, enum NNFirstNodeOptions startOption, double timeLimit, int useThreads);


//###################################################################################################################################
// EXTRA_MILEAGE
//###################################################################################################################################

Solution ExtraMileage(Instance *inst, enum EMOptions emOpt, enum EMInitType emInitType);

/* 
 * Function that applies Extra Mileage heuristic on a solution that is already not empty.
 * ARGUMENTS:
 *  sol -> solution element to which apply extra mileage
 *  nCovered -> number of elements representing the first part of the tour inside sol. The function will add all remaning nodes
 *  emOpt -> specify the way thje algorithm should operate (zero is the fastest)
 * RETURN VALUE: Number in seconds representing execution time
*/
double applyExtraMileage(Solution *sol, size_t nCovered, enum EMOptions emOpt);


//###################################################################################################################################
// 2OPT
//###################################################################################################################################

Solution _2OptBestFix(Solution *sol, enum _2OptOptions option);

double apply2OptBestFix(Solution *sol, enum _2OptOptions option);


//###################################################################################################################################
// VARIABLE_NEIGHBORHOOD
//###################################################################################################################################

Solution VariableNeighborhood(Instance *inst, enum VNSInitType config);


#endif //TSP_FUNCTIONS_H
