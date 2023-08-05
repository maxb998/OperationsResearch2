#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>

typedef struct
{
    double offset;
    int index1;
    int index2;

}SwapInformation;

static inline SwapInformation randomSwap(Solution *sol);

static inline void updateT(float *temperature);

static inline void normalizeDelta(double *deltaCost, Solution *sol);

static bool keepMove(double threshold);

void SimulatedAnnealing(Solution *sol, double timeLimit)
{
    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    // check of the input solution
    if (!checkSolution(sol))
    throwError(sol->instance, sol, "SimulatedAnnealing: Input solution is not valid");

    // initialization of variables
    int iterations = 0;
    float T = 10000;    // BETTER IF CHANGES WRT #NODES
    double offset = 0.0;
    double threshold = 0.0;

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
        iterations++;
        offset = randomSwap(newIndexPath, sol->instance);

        if (offset < 0) // then we keep the swap
        {
            for (int i = 0; i < sol->instance->nNodes; i++) sol->indexPath[i] = newIndexPath[i];
            sol->cost = computeSolutionCost(sol);
        }
        else // we implement the move with some probability
        {
            normalizeDelta(&offset, sol);
            threshold = exp(-offset/T);

            if (keepMove(threshold)) // then we keep the swap
            {
                for (int i = 0; i < sol->instance->nNodes; i++) sol->indexPath[i] = newIndexPath[i];
                sol->cost = computeSolutionCost(sol);
            }
            
            updateT(&T);
        }
    }
    
    free(newIndexPath);
}

static inline SwapInformation randomSwap(Solution *sol)
{
    int nNodes = sol->instance->nNodes;
    int *X = sol->instance->X;
    int *Y = sol->instance->Y;
    int edgeWeightType = sol->instance->params.edgeWeightType;

    int index1 = (int)((long)rand() * (long)nNodes / ((long)RAND_MAX+1L));
    int index2 = (int)((long)rand() * (long)nNodes / ((long)RAND_MAX+1L));
    {
        register int temp;
        if(index2 < index1) swapElems(index1, index2, temp);
    }// we want index1 < index2

    double oldArcsCost = 0;
    double newArcsCost = 0;
    double costEdge1 = 0;   // cost edge (index1-1, index1)
    double costEdge2 = 0;   // cost edge (index1, index1+1)
    double costEdge3 = 0;   // cost edge (index2-1, index2)
    double costEdge4 = 0;   // cost edge (index2, index2+1)

    int *path = malloc(sizeof(int) * sol->instance->nNodes);
    for (int i = 0; i < sol->instance->nNodes; i++) path[i] = sol->indexPath[i];

    costEdge2 = exactEdgeCost(X[index1], Y[index1], X[index2], Y[index2], edgeWeightType);
    
    if (index1 == 0) costEdge1 = exactEdgeCost(X[path[nNodes-1]], Y[path[nNodes-1]], X[path[index1]], Y[path[index1]], edgeWeightType);
    else costEdge1 = exactEdgeCost(X[path[index1-1]], Y[path[index1-1]], X[path[index1]], Y[path[index1]], edgeWeightType);

    if (index1 == nNodes-1) costEdge2 = exactEdgeCost(X[path[index1]], Y[path[index1]], X[path[0]], Y[path[0]], edgeWeightType);
    else costEdge2 = exactEdgeCost(X[path[index1]], Y[path[index1]], X[path[index1+1]], Y[path[index1+1]], edgeWeightType);
    
    if (index2 == 0) costEdge3 = exactEdgeCost(X[path[nNodes-1]], Y[path[nNodes-1]], X[path[index2]], Y[path[index2]], edgeWeightType);
    else costEdge3 = exactEdgeCost(X[path[index2-1]], Y[path[index2-1]], X[path[index2]], Y[path[index2]], edgeWeightType);

    if (index2 == nNodes-1) costEdge4 = exactEdgeCost(X[path[index2]], Y[path[index2]], X[path[0]], Y[path[0]], edgeWeightType);
    else costEdge4 = exactEdgeCost(X[path[index2]], Y[path[index2]], X[path[index2+1]], Y[path[index2+1]], edgeWeightType);

    oldArcsCost = costEdge1 + costEdge2 + costEdge3 + costEdge4;
    // if the nodes are adjacent we have to remove an edge that we counted twice
    if(index2 - index1 == 1 || index2 - index1 == nNodes-1) oldArcsCost -= exactEdgeCost(X[path[index1]], Y[path[index1]], X[path[index2]], Y[path[index2]], edgeWeightType);

    {
        register int temp;
        swapElems(path[index1], path[index2], temp);
    }

    if (index1 == 0) costEdge1 = exactEdgeCost(X[path[nNodes-1]], Y[path[nNodes-1]], X[path[index1]], Y[path[index1]], edgeWeightType);
    else costEdge1 = exactEdgeCost(X[path[index1-1]], Y[path[index1-1]], X[path[index1]], Y[path[index1]], edgeWeightType);

    if (index1 == nNodes-1) costEdge2 = exactEdgeCost(X[path[index1]], Y[path[index1]], X[path[0]], Y[path[0]], edgeWeightType);
    else costEdge2 = exactEdgeCost(X[path[index1]], Y[path[index1]], X[path[index1+1]], Y[path[index1+1]], edgeWeightType);
    
    if (index2 == 0) costEdge3 = exactEdgeCost(X[path[nNodes-1]], Y[path[nNodes-1]], X[path[index2]], Y[path[index2]], edgeWeightType);
    else costEdge3 = exactEdgeCost(X[path[index2-1]], Y[path[index2-1]], X[path[index2]], Y[path[index2]], edgeWeightType);

    if (index2 == nNodes-1) costEdge4 = exactEdgeCost(X[path[index2]], Y[path[index2]], X[path[0]], Y[path[0]], edgeWeightType);
    else costEdge4 = exactEdgeCost(X[path[index2]], Y[path[index2]], X[path[index2+1]], Y[path[index2+1]], edgeWeightType);

    newArcsCost = costEdge1 + costEdge2 + costEdge3 + costEdge4;
    // if the nodes are adjacent we have to remove an edge that we counted twice
    if(index2 - index1 == 1 || index2 - index1 == nNodes-1) newArcsCost -= exactEdgeCost(X[path[index1]], Y[path[index1]], X[path[index2]], Y[path[index2]], edgeWeightType);

    SwapInformation swapInfo = {.offset = newArcsCost - oldArcsCost, .index1 = index1, .index2 = index2};
    free(path);

    return swapInfo;
}

static inline void updateT(float *temperature)
{
    *temperature = 0.99 * (*temperature);
}

static inline void normalizeDelta(double *offset, Solution *sol)
{
    *offset = *offset / (sol->cost / sol->instance->nNodes);
}

static inline bool keepMove(double threshold)
{
    if(threshold*RAND_MAX > rand()) return true;
    else return false;
}