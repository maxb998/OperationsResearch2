#ifndef TSP_H
#define TSP_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>

#include <pthread.h>    // for multithreadin
#include <immintrin.h>  // for avx simd instructions

#define SMALLX 1e-6
#define EPSILON 1e-9

// data structures
typedef struct
{
    //int computeDistFirst;
    int randomSeed;
    char inputFile[1000];
    int threadsCount;
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
    //double *xNodes;
    //double *yNodes;
    double *coords;     // all x first and then the y
    double *edgeCost;   // matrix with the cost of all edges (can use -1 if edge does not exists)

    Parameters params;
    
    GlobalData global;
    
} Instance;



#endif //TSP_H