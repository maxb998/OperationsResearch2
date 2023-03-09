#ifndef TSP_H
#define TSP_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>

#include <pthread.h>    // for multithreadin
#include <immintrin.h>  // for avx simd instructions

//#define VERBOSE 10
const int VERBOSE = 150; // printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)
const double smallX = 1e-6;
const double epsilon = 1e-9;

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