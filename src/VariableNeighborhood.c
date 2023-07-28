#include "Tsp.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time
#include <stdio.h>


// Returns 5 random different integers in [0,nNodes) making sure they differ of more than 1 unit
static inline int * random5Edges(int nNodes, Solution * sol);

static inline void _5Kick(Solution *sol);

static inline void sort5int(int * array);

void VariableNeighborhood(Solution *sol, double timeLimit, int nThreads)
{
    // time limit management
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    if ((nThreads < 0) || (nThreads > MAX_THREADS))
        throwError(sol->instance, sol, "VariableNeighborhood: nThreads value is not valid: %d", nThreads);
    else if (nThreads == 0)
        nThreads = sol->instance->params.nThreads;

    Instance *inst = sol->instance;

    if(inst->nNodes <= 9)
        throwError(inst, sol, "VariableNeighborhood: number of nodes too low to run VNS. #Nodes: ", inst->nNodes);

    if (!checkSolution(sol))
        throwError(inst, sol, "VariableNeighborhood: Input solution is not valid");

    apply2OptBestFix(sol, _2OPT_AVX_ST);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);

    while (currentTime < startTime + timeLimit)
    {
        _5Kick(sol);
        apply2OptBestFix(sol, _2OPT_AVX_ST);

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
        currentTime = cvtTimespec2Double(timeStruct);
    }
}

static inline void _5Kick(Solution *sol)
{

    // Contains the index to 5 random non consecutive edges, sorted from smaller index to larger
    int * randomEdges = random5Edges((int)sol->instance->nNodes, sol);
    sort5int(randomEdges);

    int * newPath = malloc((sol->instance->nNodes+1)*sizeof(int));
    int index = 0;

    newPath[index] = sol->indexPath[randomEdges[0]];
    index++;
    for(int i = randomEdges[1]+1; i <= randomEdges[2]; i++)
    {
        newPath[index] = sol->indexPath[i];
        index++;
    }
    for(int i = randomEdges[3]+1; i <= randomEdges[4]; i++)
    {
        newPath[index] = sol->indexPath[i];
        index++;
    }
    for(int i = randomEdges[0]+1; i <= randomEdges[1]; i++)
    {
        newPath[index] = sol->indexPath[i];
        index++;
    }
    for(int i = randomEdges[2]+1; i <= randomEdges[3]; i++)
    {
        newPath[index] = sol->indexPath[i];
        index++;
    }
    for(int i = randomEdges[4]+1; i <= sol->instance->nNodes; i++)
    {
        newPath[index] = sol->indexPath[i];
        index++;
    }
    for(int i = 1; i <= randomEdges[0]; i++)
    {
        newPath[index] = sol->indexPath[i];
        index++;
    }
    
    LOG(LOG_LVL_DEBUG, "_5Kick -> Edges mixed: [%d-%d], [%d-%d], [%d-%d], [%d-%d], [%d-%d]", sol->indexPath[randomEdges[0]], sol->indexPath[randomEdges[0]+1], sol->indexPath[randomEdges[1]], sol->indexPath[randomEdges[1]+1], sol->indexPath[randomEdges[2]], sol->indexPath[randomEdges[2]+1], sol->indexPath[randomEdges[3]], sol->indexPath[randomEdges[3]+1], sol->indexPath[randomEdges[4]], sol->indexPath[randomEdges[4]+1]);

    for(int i = 0; i <= sol->instance->nNodes; i++)
    {
        sol->indexPath[i] = newPath[i];
        sol->X[i] = sol->instance->X[newPath[i]];
        sol->Y[i] = sol->instance->Y[newPath[i]];
    }
    
    sol->cost = computeSolutionCost(sol);

    free(randomEdges);
    free(newPath);
}

static inline int * random5Edges(int nNodes, Solution * sol)
{
    int * randomEdges = malloc(5*sizeof(int));
    for(int i = 0; i < 5; i++)
    {
        randomEdges[i] = rand()%nNodes;
        for(int j = 0; j < i; j++)
        {
            if((randomEdges[i] == randomEdges[j]) || (abs(randomEdges[i] - randomEdges[j]) == 1) || (abs(randomEdges[i] - randomEdges[j]) == sol->instance->nNodes-1)) 
            {
                i--;
                break;            
            }
        }
    }
    return randomEdges;
}

static inline void sort5int(int * array)
{	
	for(int j = 0; j < 4; j++)
	{
		int min = j;
		for(int i = j+1; i < 5; i++)
			if(array[i] < array[min])
                min = i;
		register int temp;
        swapElems(array[j], array[min], temp);
	}
}