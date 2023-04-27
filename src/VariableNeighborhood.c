#include "VariableNeighborhood.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"
#include "NearestNeighbor.h"
#include "ExtraMileage.h"
#include "2Opt.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

/* Checks wheter the time passed since the initialization of start has passed timeLimit.
    Returns 1 in that case, 0 otherwise.*/
static inline int checkTime(struct timespec start, double timeLimit);

// Returns 5 random different integers in [0,nNodes) making sure they differ of more than 1 unit
static inline int * random5Edges(int nNodes, Solution * sol);

static inline void _5Kick(Solution *sol);

static inline void sort5int(int * array);

Solution VariableNeighborhood(Instance *inst, enum VNSInitType config)
{
    // time limit management
    double placeholderTime = 10000.0;
    struct timespec start;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    Solution sol = newSolution(inst);

    // check if a seed for random has been passed as argument
    if(inst->params.randomSeed != -1) srand(inst->params.randomSeed);
    else throwError(inst, &sol, "VariableNeighborhood: random seed has not been passed as argument");
    
    if((config != VNS_INIT_NN) && (config != VNS_INIT_EM)) throwError(inst, &sol, "VariableNeighborhood: incorrect argument for config");
    else 
    {
        if(config == VNS_INIT_NN) // Find local minimum with Nearest Neighbour
        {
            // Compute a solution with Nearest Neighbour and optimize it with 2-opt
            sol = NearestNeighbor(inst, inst->params.nnFirstNodeOption, inst->params.tlim/20, 1);
            apply2OptBestFix(&sol, _2OPT_AVX_ST);
            while(checkTime(start, placeholderTime) == 0)
            {
                _5Kick(&sol);
                apply2OptBestFix(&sol, _2OPT_AVX_ST);
            }

        }else // (config == VNS_INIT_EM)
        {   
            // Compute a solution with Extra Mileage and optimize it with 2-opt
            sol = ExtraMileage(inst, EM_OPTION_AVX, EM_INIT_RANDOM);
            apply2OptBestFix(&sol, _2OPT_AVX_ST);
            while(checkTime(start, placeholderTime) == 0)
            {
                _5Kick(&sol);
                apply2OptBestFix(&sol, _2OPT_AVX_ST);
                
            }
            
        }
    }
    return sol;
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

    /*
    for (size_t i = 0; i <= sol->instance->nNodes; i++)
    {
        printf("%d: %d -> %d\n", i+1, sol->indexPath[i], newPath[i]);
    }
    */
    
    LOG(LOG_LVL_DEBUG, "_5Kick -> Edges mixed: [%d-%d], [%d-%d], [%d-%d], [%d-%d], [%d-%d]", sol->indexPath[randomEdges[0]], sol->indexPath[randomEdges[0]+1], sol->indexPath[randomEdges[1]], sol->indexPath[randomEdges[1]+1], sol->indexPath[randomEdges[2]], sol->indexPath[randomEdges[2]+1], sol->indexPath[randomEdges[3]], sol->indexPath[randomEdges[3]+1], sol->indexPath[randomEdges[4]], sol->indexPath[randomEdges[4]+1]);

    for(size_t i = 0; i <= sol->instance->nNodes; i++)
    {
        sol->indexPath[i] = newPath[i];
        sol->X[i] = sol->instance->X[newPath[i]];
        sol->Y[i] = sol->instance->Y[newPath[i]];
    }

    free(randomEdges);
    free(newPath);
}

static inline int checkTime(struct timespec start, double timeLimit)
{
    struct timespec currentTime;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currentTime);
    double elapsed = ((currentTime.tv_sec - start.tv_sec) + (currentTime.tv_nsec - start.tv_nsec) / 1000000000.0); // nsed are useless if we are counting time in minutes
    
    if(elapsed < timeLimit) return 0;
    else return 1;
}

static inline int * random5Edges(int nNodes, Solution * sol)
{
    if(nNodes <= 9) throwError(sol->instance, sol, "random5edges: number of nodes too low to run VNS. #Nodes: ", nNodes);
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
    int temp;	
	for(int j = 0; j < 4; j++)
	{
		int min = j;
		for(int i = j+1; i < 5; i++)
		{
			if(array[i] < array[min]) min = i;
		}
		temp = array[j];
		array[j] = array[min];
		array[min] = temp;
	}
}