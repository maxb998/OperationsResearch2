#ifndef TSP_INCLUDES
#define TSP_INCLUDES

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <stdint.h>

#include <stdarg.h> // used for logger va_list
#include <getopt.h> // args parsing by POSIX
#include <limits.h> // used to get hard limits mainly to do error checking(like overflows and such)
#include <time.h>	// to manage the time limits for the meta-heuristics

#include <pthread.h>    // for multithreadin
#include <immintrin.h>  // for avx simd instructions

#endif //TSP_INCLUDES


#ifndef TSP_DATA_BASE
#define TSP_DATA_BASE

#define MAX_THREADS 32
#define USE_APPROXIMATED_DISTANCES 1
#define LOG_LEVEL LOG_LVL_EVERYTHING

// size of avx vector. 4 is vector of doubles 64bits, 8 is vector of floats 32bits
#define AVX_VEC_SIZE 8

#define SMALLX 1e-6
#define EPSILON 1e-9

enum edgeWeightType{
	EUC_2D, // euclidean distance 2d
	MAN_2D, // manhattan distance 2d
	MAX_2D, // maximum distance 2d
	CEIL_2D, // euclidean 2d rounded up
	ATT, // special distance for problems att48 and att532
	EUC_3D, // euclidean distance 3d
	MAN_3D, // manhattan distance 3d
	MAX_3D, // maximum distance 3d
    GEO, // geographical distance
	XRAY1, // special distance for crystallography problems v1
	XRAY2, // special distance for crystallography problems v2
	EXPLICIT, // weights are specified in the file
	SPECIAL // special type of distance documented elsewhere
};

// data structures
typedef struct
{
    int edgeWeightType;
    int randomSeed;			// if no value has been specified as argument GRASP won't be enabled
	int roundWeights;
    char inputFile[1000];
	char name[200];
    int nThreads;	// if no value has been specified as argument its default value is the number of processors in the machine
} Parameters;

typedef struct
{
    // data
    size_t nNodes;
    float *X;
    float *Y;

    Parameters params;
	void *bestSol;
} Instance;

typedef struct
{
    float bestCost;    // best solution found cost
	int bestCostRounded;

	double execTime;

	// the solution
	float *X;
	float *Y;

	Instance *instance;
} Solution;

#endif //TSP_DATA_BASE

#ifndef TSP_UTILITIES
#define TSP_UTILITIES

#define GRASP_COEFF (int)RAND_MAX*0.9

enum logLevel{
	LOG_LVL_ERROR, // 0
	LOG_LVL_CRITICAL, // 1
	LOG_LVL_WARNING, // 2
	LOG_LVL_NOTICE, // 3
	LOG_LVL_LOG, // 4
	LOG_LVL_DEBUG, // 5
	LOG_LVL_EVERYTHING // 6
};

// Returns initialized/empty instance
Instance newInstance ();

// Returns initialized solution struct for the specified Instance pointed by inst
Solution newSolution (Instance *inst);

// Frees the allocated memory in the Instance pointed by d
void destroyInstance (Instance *inst);

// Frees the allocated memory in the Solution pointed by s
void destroySolution (Solution *sol);

/*Function used to log informations with level.
* lvl -> logLevel: Desired level of logging priority. If this value is greater than global value than LOG does not print anything
* line, ... -> string, args: Arguments that are directly fed to printf
*/
int LOG (enum logLevel lvl, char * line, ...);

// Launch fatal error, free memory and exit with code 1
void throwError (Instance *inst, Solution *sol, char * line, ...);

// Parse commandline arguments stored in argv and save relevant information to d->params
void parseArgs (Instance *inst, int argc, char *argv[]);

/* Checks that:
	1. First and last node in the path are the same (closed circuit)
	2. Nodes are not getting covered more than once
	3. All nodes are covered along the path*/
int solutionCheck(Instance *inst);

// Check the correctness of the cost of the solution stored in Instance inst.
int costCheck(Instance *inst);

// Read file with .tsp extension according to tsplib specifications, complete with file sintax error checking
void readFile (Instance *inst);

void saveSolution(Solution *sol);

/*Plot solution using gnuplot. Does NOT check for errors on input
* d	-> Instance to plot
* plotPixelSize	-> string: Plot window size in pixel specified with format: "<WIDTH>,<HEIGHT>"
* pointColor -> string: Color of the circle representing the point, eg "black" or "red"
* tourPointColor -> string: Color of the 'X' on top of the point circle of color pointColor. Format and types is the same as for pointColor
* pointSize -> int: Size of the points
*/
void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize);

double computeSquaredCost_VEC(Solution *sol);

double computeSquaredCost(Solution *sol);

#endif //TSP_UTILITIES


#ifndef NEAREST_NEIGHBOUR
#define NEAREST_NEIGHBOUR

// Computes the Nearest Neighbour heuristic starting from every node, and saving the path with minimum cost into the Instance
// Creates a thread for every logic processor in the machine
Solution NearestNeighbour(Instance *inst);

#endif //NEAREST_NEIGHBOUR


#ifndef EXTRA_MILEAGE
#define EXTRA_MILEAGE

double solveExtraMileage(Instance *inst);

#endif //EXTRA_MILEAGE


#ifndef _2OPT
#define _2OPT

double _2optBestFix(Instance *inst);

#endif //_2OPT

#ifndef VARIABLE_NEIGHBORHOOD
#define VARIABLE_NEIGHBORHOOD

double VariableNeighborhood(Instance *inst, int configuration);

#endif // VARIABLE_NEIGHBORHOOD

