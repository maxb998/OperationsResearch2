#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

static inline double randomSwap(int *path, Solution *sol);

static inline void updateT(float *temperature);

static inline void normalizeDelta(double *deltaCost, Solution *sol);

static bool keepMove(double *threshold);

void SimulatedAnnealing(Solution *sol, double timeLimit)
{
    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // check of the input solution
    if (!checkSolution(sol))
    throwError(sol->instance, sol, "SimulatedAnnealing: Input solution is not valid");

    double currentCost = sol->cost;
    int iterations = 0;
    float T = 10000;    // BETTER IF CHANGES WRT #NODES
    int pathSize = sol->instance->nNodes;

    // Contains a copy of the current indexPath, which will be modified and from which we will update the solution
    int *newIndexPath = malloc(sizeof(int) * sol->instance->nNodes);
    for (int i = 0; i < sol->instance->nNodes; i++)
    {
        newIndexPath[i] = sol->indexPath[i];
    }

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    while (currentTime < startTime + timeLimit)
    {

    }
    
    free(newIndexPath);
}

static inline double randomSwap(int *path, Solution *sol)
{
    int nNodes = sol->instance->nNodes;
    int index1 = (int)((long)rand() * (long)nNodes / ((long)RAND_MAX+1L));
    int index2 = (int)((long)rand() * (long)nNodes / ((long)RAND_MAX+1L));
    {
        register int temp;
        if(index2 < index1) swapElems(index1, index2, temp);
    }// we want index1 < index2

    double costEdge1 = 0;   // cost edge (index1-1, index1)
    double costEdge2 = 0;   // cost edge (index1, index1+1) or (index1,index2) if the nodes are adjacent
    double costEdge3 = 0;   // cost edge (index2-1, index2) or 0 if the nodes are adjacent
    double costEdge4 = 0;   // cost edge (index2, index2+1)

    
    if(index2-index1 == 1 || index2-index1 == nNodes-1) // nodes are adjacent
    {
        costEdge2 = exactEdgeCost(sol->X[index1], sol->Y[index1], sol->X[index2], sol->Y[index2], sol->instance->params.edgeWeightType);
    }else
    {
        costEdge2 = exactEdgeCost(sol->X[index1], sol->Y[index1], sol->X[index1+1], sol->Y[index1+1], sol->instance->params.edgeWeightType);
        costEdge3 = exactEdgeCost(sol->X[index2-1], sol->Y[index2-1], sol->X[index2], sol->Y[index2], sol->instance->params.edgeWeightType);
    }

    if(index1 == 0) costEdge1 = exactEdgeCost(sol->X[nNodes-1], sol->Y[nNodes-1], sol->X[index1], sol->Y[index1], sol->instance->params.edgeWeightType);
    else costEdge1 = exactEdgeCost(sol->X[index1-1], sol->Y[index1-1], sol->X[index1], sol->Y[index1], sol->instance->params.edgeWeightType);


        
    {
        register int temp;
        swapElems(path[index1], path[index2], temp);
    }

    double costOffset = 0;

    return costOffset;
}

static inline void updateT(float *temperature)
{
    *temperature = 0.99 * (*temperature);
}

static inline void normalizeDelta(double *deltaCost, Solution *sol)
{
    *deltaCost = *deltaCost / (sol->cost / sol->instance->nNodes);
}

static bool keepMove(double *threshold)
{
    if(*threshold*RAND_MAX > rand()) return true;
    else return false;
}