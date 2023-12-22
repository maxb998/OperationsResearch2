#ifndef TSP_FUNCTIONS_H
#define TSP_FUNCTIONS_H

#include "TspBase.h"
#include "EdgeCostFunctions.h"

#include "limits.h"

//#define DEBUG

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
void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize, const bool printIndex);

//###################################################################################################################################
// UTILITIES
//###################################################################################################################################

// Swap elem1 and elem2. 
#define swapElems(elem1,elem2) { register typeof(elem1) swapVarTemp = elem1; elem1 = elem2; elem2 = swapVarTemp; }

// Convert timespec struct to a double time in seconds
#define cvtTimespec2Double(t) (double)t.tv_sec + (double)t.tv_nsec / 1000000000.0

// Use rand_r to generate a random integer in the range [min,max) while avoiding low entropy bits (statiscally better than doing "rand() % (n+1)") https://codereview.stackexchange.com/questions/159604/uniform-random-numbers-in-an-integer-interval-in-c
#define genRandom(rndStatePtr,min,max) (int)((long)rand_r(rndStatePtr) * (long)((max)-(min)) / (RAND_MAX + 1L)) + (min)

// Convert float to fixed point 128 bit
#define cvtFloat2Cost(fpVar) ((__uint128_t)((double)fpVar * (double)ULONG_MAX))

// Convert fixed point 128bit to double
#define cvtCost2Double(cost) ((double)cost / (double)ULONG_MAX)

/*!
* @brief Generate empty instance
* @result Empty Instance
*/
Instance newInstance ();

/*!
* @brief Initialize a solution struct for the specified instance
* @param inst Instance related to the new solution
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
* @brief Fully clones the passed Solution. Pointers are not copied over, but every element inside indexPath is copied from as well as all other data.
* @param src Solution to clone
* @param dst Destination solution (the size of the instance must be the same)
*/
void cloneSolution(Solution *src, Solution *dst);

/*!
* @brief Set the log level value which indicates how many information shall be printed as output on the console.
* @param lvl set log level desired
*/
void setLogLevel(enum LogLevel lvl);

/*!
* @brief Print a message on output
* @param lvl Desired level of logging priority. If this value is greater than global value than LOG does not print anything.
* @param line Message with by printf format args
*/
void LOG (enum LogLevel lvl, char * line, ...);

/*!
* @brief Print a message with LOG_LVL_ERROR and then exit the program.
* @param line Error message with by printf format args
*/
void throwError (char * line, ...);

/*!
* @brief Checks if a solution is feasible both in the sequences of sol.{x,y,indexPath} and in the value of the cost.
* @param sol Solution to check.
* @returns true if solution is correct, false otherwise
*/
bool checkSolution(Solution *sol);

/*!
* @brief Compute/Recompite the cost of the solution given.
* @param sol Solution of which compute the cost.
* @result cost of solution
*/
__uint128_t computeSolutionCost(Solution *sol);

/*!
* @brief Perform argsort using quicksort algorithm with random pivot
* @param arr Floating point array to sort.
* @param indexes Array where to store the sorted indexes of arr. (memory MUST be preallocated (n*sizeof(int)))
* @param n Amount of elements in arr and indexes
*/
void argsort(float *arr, int *indexes, int n);

/*!
* @brief Sort arr using quicksort algorithm with random pivot
* @param arr Floating point array to sort.
* @param n Amount of elements in arr
*/
void sort(float *arr, int n);

#if (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
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

#endif

//###################################################################################################################################
// NEAREST_NEIGHBOR
//###################################################################################################################################

/*!
* @brief Find a the best Solution using Nearest Neighbor Heuristic
* @param inst Instance to solve.
* @param timeLimit Time Limit.
* @result Solution obtained using Nearest Neighbor.
*/
Solution NearestNeighbor(Instance *inst, double timeLimit);


//###################################################################################################################################
// EXTRA_MILEAGE
//###################################################################################################################################

/*!
* @brief Find a the best Solution using Extra Mileage Heuristic.
* @param inst Instance to solve.
* @param timeLimit Time limit.
* @result Solution obtained using Extra Mileage.
*/
Solution ExtraMileage(Instance *inst, double timeLimit);

//###################################################################################################################################
// 2OPT
//###################################################################################################################################
/*!
* @brief Set the option to view performance related statistics on 2opt runs
# @param val True to enable, false to disable
*/
void set2OptPerformanceBenchmarkLog(bool val);

/*!
* @brief  Applies 2Opt solution optimizer to sol.
* @param sol Solution to optimize.
*/
void apply2OptBestFix(Solution *sol);


#if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
/*!
* @brief  Same as apply2OptBestFix, but expects costCache array and, if using AVX, X and Y arrays all coherent with sol.indexPath. If not using AVX pass X = Y = NULL
* @param sol Solution to optimize.
*/
void apply2OptBestFix_fastIteratively(Solution *sol, float *X, float *Y, float *costCache);
#elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
/*!
* @brief  Same as apply2OptBestFix, but expects costCache array and, if using AVX, X and Y arrays all coherent with sol.indexPath. If not using AVX pass X = Y = NULL
* @param sol Solution to optimize.
*/
void apply2OptBestFix_fastIteratively(Solution *sol, float *costCache);
#endif

//###################################################################################################################################
// TABU_SEARCH
//###################################################################################################################################

/*!
* @brief Runs Tabu Search on input solution sol.
* @param sol Starting point solution.
* @param timeLimit Time limit.
*/
void TabuSearch(Solution *sol, double timeLimit);

//###################################################################################################################################
// VARIABLE_NEIGHBORHOOD
//###################################################################################################################################

/*!
* @brief  Run Variable Neighborhood Search on inst and return the best solution found.
* @param sol Starting point solution for the VNS algorithm.
* @param timeLimit Time limit.
*/
void VariableNeighborhoodSearch(Solution *sol, double timeLimit);


//###################################################################################################################################
// SIMULATED_ANNEALING
//###################################################################################################################################

/*!
* @brief Runs Simulated Annealing on input solution sol.
* @param sol Starting point solution.
* @param timeLimit Time limit.
*/
void SimulatedAnnealing(Solution *sol, double timeLimit);



//###################################################################################################################################
// GENETIC_ALGORITHM
//###################################################################################################################################

/*!
* @brief Runs Genetic Algorithm with population, crossover, population and reintroduction amount specified in inst.params
* @param inst Instance
* @param timeLimit Time limit.
* @result Best solution inside the population after the time limit.
*/
Solution GeneticAlgorithm(Instance *inst, double timeLimit);

#endif //TSP_FUNCTIONS_H
