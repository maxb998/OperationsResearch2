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

#define SMALLX 1e-6
#define EPSILON 1e-9

enum edgeWeightType{
	EUC_2D, // 0, euclidean distance 2d
	MAN_2D, // 1, manhattan distance 2d
	MAX_2D, // 2, maximum distance 2d
	CEIL_2D, // 3, euclidean 2d rounded up
    GEO, // 4,
	ATT, // 5, special distance for problems att48 and att532
	XRAY1, // 6, special distance for crystallography problems v1
	XRAY2, // 7, special distance for crystallography problems v2
	EUC_3D, // 8, euclidean distance 3d
	MAN_3D, // 9, manhattan distance 3d
	MAX_3D, // 10, maximum distance 3d
	EXPLICIT, // 11, weights are specified in the file
	SPECIAL // 12, special type of distance documented elsewhere
};

// data structures
typedef struct
{
    //int computeDistFirst;
    int edgeWeightType;
    int randomSeed;
    char inputFile[1000];
    size_t threadsCount;
    int useAVX;
} Parameters;

typedef struct
{
    double bestCost;    // best solution found cost
    int *bestSolution;  // array containing sequence of nodes representing the optimal solution
} GlobalData;


typedef struct
{
    // data
    size_t nodesCount;
    double *X;
    double *Y;
    //double *coords;     // all x first and then the y
    double *edgeCost;   // matrix with the cost of all edges (can use -1 if edge does not exists)

    Parameters params;
    
    GlobalData global;
    
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

int LOG (enum logLevel lvl, char * line, ...);

void parseArgs (Instance *d, int argc, char *argv[]);

void readFile (Instance *d);

#endif //TSP_UTILITIES

#ifndef TSP
#define TSP

#define MAX_THREADS 255
#define AVX_VEC_SIZE 4

int computeSquaredDistanceMatrix(Instance *d);

#endif //TSP