#ifndef TSP_INCLUDES
#define TSP_INCLUDES

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include <stdarg.h> // used for logger va_list
#include <getopt.h> // args parsing by POSIX
#include <limits.h> // used to get hard limits mainly to do error checking(like overflows and such)

#include <pthread.h>    // for multithreadin
#include <immintrin.h>  // for avx simd instructions

#endif //TSP_INCLUDES


#ifndef TSP_DATA_STRUCTURES
#define TSP_DATA_STRUCTURES

#define MAX_THREADS 255

// size of avx vector. 4 is vector of doubles 64bits, 8 is vector of floats 32bits
#define AVX_VEC_SIZE 8

#define DOUBLE_MAX 1.79769313486231570e+308
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
    int randomSeed;
	int roundWeights;
    char inputFile[1000];
	char name[200];
    size_t threadsCount;
} Parameters;

typedef struct
{
    float bestCost;    // best solution found cost
	int bestRoundedSol;
    int *bestSolution;  // array containing sequence of nodes representing the optimal solution
} Solution;

typedef struct
{
	// if we require rounded weights mat will be NULL and roundedMat will point to the allocated matrix
	float * mat;
	int * roundedMat;
	size_t rowSizeMem;
} EdgeCostMatStruct;


typedef struct
{
    // data
    size_t nodesCount;
    float *X;
    float *Y;
    //double *coords;     // all x first and then the y
    EdgeCostMatStruct edgeCost;   // matrix with the cost of all edges (can use -1 if edge does not exists)

    Parameters params;
    
    Solution solution;
    
} Instance;

#endif //TSP_DATA_STRUCTURES

#ifndef TSP_UTILITIES
#define TSP_UTILITIES

enum logLevel{
	LOG_LVL_ERROR, // 0
	LOG_LVL_CRITICAL, // 1
	LOG_LVL_WARNING, // 2
	LOG_LVL_NOTICE, // 3
	LOG_LVL_LOG, // 4
	LOG_LVL_DEBUG, // 5
	LOG_LVL_EVERYTHING // 6
};

void initInstance(Instance *d);
void freeInstance(Instance *d);

int LOG (enum logLevel lvl, char * line, ...);

// launch fatal error, free memory and exit with code 1
void throwError (Instance *d, char * line, ...);

void parseArgs (Instance *d, int argc, char *argv[]);

void readFile (Instance *d);

void saveSolution(Instance *d);

void plotSolution(Instance *d);

#endif //TSP_UTILITIES

#ifndef DISTANCE_MATRIX
#define DISTANCE_MATRIX

void printDistanceMatrix(Instance *d, int showEndRowPlaceholder);

int computeSquaredDistanceMatrix(Instance *d);

#endif //DISTANCE_MATRIX