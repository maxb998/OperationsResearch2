#ifndef TSP_FUNCTIONS_H
#define TSP_FUNCTIONS_H

#include "TspBase.h"
#include "EdgeCostFunctions.h"

//###################################################################################################################################
// ARG_PARSER 
//###################################################################################################################################

/*!
* @brief Parse arguments and store them in inst->params
* @param inst Pointer initialized instance(newInstance())
*/
void argParse(Instance * inst, int argc, char *argv[]);

/*!
* @brief Print parameters and defaults value for the run. Call this after argParse
* @param inst Instance
*/
void printInfo(Instance *inst);

//###################################################################################################################################
// TSP_IOUTILS
//###################################################################################################################################

/*!
* @brief Read .tsp file according to tsplib specification. Very strict formatting must be applied to .tsp file
* @param inst Instance to which save data
* @result Time to read file
*/
double readFile (Instance *inst);

/*!
* @brief Save solution inside run/ directory with all parameters used
* @param inst Instance to which save data
* @result Time to read file
*/
void saveSolution(Solution *sol, int argc, char *argv[]);

/*!
* @brief Plot solution using gnuplot. Does NOT check for errors on input
* @param sol Solution to plot
* @param plotPixelSize Plot window size in pixel specified with format: "<WIDTH>,<HEIGHT>"
* @param pointColor Color of the circle representing the point, eg "black" or "red"
* @param tourPointColor Color of the 'X' on top of the point circle of color pointColor. Format and types is the same as for pointColor
* @param pointSize Size of the points
* @param printIndex Set to 1 to print index of each point as label on plot
*/
void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize, const int printIndex);

//###################################################################################################################################
// UTILITIES
//###################################################################################################################################

// Swap elem1 and elem2. Can be done with any type of variable, however a temporary variable "tmp" of the same type of elem1 and elem2 MUST be provided.
#define swapElems(elem1,elem2,tmp) tmp = elem1; elem1 = elem2; elem2 = tmp

/*!
* @brief Generate empty instance
* @result Empty Instance
*/
Instance newInstance ();

/*!
* @brief Initialize a solution struct for the specified instance
* @param inst 
* @result Empty Solution with allocated memory
*/
Solution newSolution (Instance *inst);

/*!
* @brief Frees the allocated memory in the passed instance
* @param inst Instance which memory needs to be freed
*/
void destroyInstance (Instance *inst);

/*!
* @brief Frees the allocated memory in the passed solution (not on the instance pointed inside the solution)
* @param sol Solution which memory needs to be freed
*/
void destroySolution (Solution *sol);

/*!
* @brief Fully clones the passed Solution. Pointers are not just copied over, they point to new memory allocation that will contain a copy of the memory pointed by the pointers in sol.
* @param sol Solution to clone
* @result Cloned Solution
*/
Solution cloneSolution(Solution *sol);

/*!
* @brief Set the log level value which indicates how many information shall be printed as output on the console.
* @param lvl set log level desired
*/
void setLogLevel(enum logLevel lvl);

/*!
* @brief Print a message on output
* @param lvl Desired level of logging priority. If this value is greater than global value than LOG does not print anything.
* @param line Arguments that are directly fed to printf
*/
void LOG (enum logLevel lvl, char * line, ...);

/*!
* @brief Print a message with LOG_LVL_ERROR and then destroys the passed Instance and Solution when available.
* @param inst Pointer to the Instance to destroy. Can be NULL.
* @param sol Pointer to the Solution to destroy. Can be NULL.
* @param line Error message.
*/
void throwError (Instance *inst, Solution *sol, char * line, ...);

/*!
* @brief Checks if a solution is feasible both in the sequences of sol.{x,y,indexPath} and in the value of the cost.
* @param sol Solution to check.
* @throw calls throwError(sol, sol.instance, ...) in case of failure
*/
void checkSolution(Solution *sol);

/*!
* @brief Compute/Recompite the cost of the solution given.
* @param sol Solution of which compute the cost.
* @result cost of solution
*/
double computeSolutionCost(Solution *sol);

/*!
* @brief Compute/Recompite the cost of the solution given using vectorization (Didn't test if this is any faster than the other one)
* @param sol Solution of which compute the cost.
* @result cost of solution
*/
double computeSolutionCostVectorized(Solution *sol);


//###################################################################################################################################
// COST_MATRIX
//###################################################################################################################################

/*!
* @brief Allocate memory for the cost matrix and compute the cost matrix.
* @param inst Instance of which compute the cost matrix.
* @result Time needed to compute the cost matrix.
*/
double computeCostMatrix(Instance *inst);

/*!
* @brief Print the cost matrix of an instance in console.
* @param inst Instance containing the pointer for the cost matrix.
* @attention Only use with very small Instances since the output is quadratic to the number of nodes
*/
void printCostMatrix(Instance *inst);

//###################################################################################################################################
// NEAREST_NEIGHBOR
//###################################################################################################################################

/*!
* @brief Find a the best Solution using Nearest Neighbor Heuristic withing a specified time limit.
* @param inst Instance to solve.
* @param startOption Defines the way in which the first point has to be chosen
* @param tlim Time Limit
* @param useThread Set to 1 to use multithreading. Set to 0 for single thread.
* @result Solution obtained using Nearest Neighbor.
*/
Solution NearestNeighbor(Instance *inst, enum NNFirstNodeOptions startOption, double tlim, int useThreads);


//###################################################################################################################################
// EXTRA_MILEAGE
//###################################################################################################################################

/*!
* @brief Find a the best Solution using Extra Mileage Heuristic withing a specified time limit.
* @param inst Instance to solve.
* @param emOpt Defines the way in which the first point has to be chosen.
* @param emInitType Defines the way in which the Solution is initialized.
* @param tlim Time limit.
* @result Solution obtained using Extra Mileage.
*/
Solution ExtraMileage(Instance *inst, enum EMOptions emOpt, enum EMInitType emInitType);

/*!
* @brief Applies Extra Mileage heuristic on a solution that is already not empty.
* @param sol Solution with first batch already inserted up to nCovered.
* @param nCovered Number of elements representing the first part of the tour inside sol. The function will add all remaning nodes.
* @param emOpt Defines the way in which the first point has to be chosen.
* @result Number in seconds representing execution time
*/
double applyExtraMileage(Solution *sol, size_t nCovered, enum EMOptions emOpt);


//###################################################################################################################################
// 2OPT
//###################################################################################################################################

/*!
* @brief Applies 2Opt solution optimizer to a cloned solution which is returned as output.
* @param sol Solution to optimize.
* @param option Type of algorithm to use. Only affects speed of computation. To use the fastest one set to 0.
* @result 2Optimized version of sol.
*/
Solution _2OptBestFix(Solution *sol, enum _2OptOptions option);

/*!
* @brief  Applies 2Opt solution optimizer to sol.
* @param sol Solution to optimize.
* @param option Type of algorithm to use. Only affects speed of computation. To use the fastest one set to 0.
* @result Number in seconds representing execution time
*/
double apply2OptBestFix(Solution *sol, enum _2OptOptions option);


//###################################################################################################################################
// VARIABLE_NEIGHBORHOOD
//###################################################################################################################################

/*!
* @brief  Run Variable Neighborhood Search on inst and return the best solution found.
* @param inst Instance to solve.
* @param config Type of Solver to use for the first solution. Choose between Nearest Neighbor and Extra Mileage.
* @result Best solution found with VNS.
*/
Solution VariableNeighborhood(Instance *inst, enum VNSInitType config);


#endif //TSP_FUNCTIONS_H