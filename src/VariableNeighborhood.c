#include "VariableNeighborhood.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"
#include "NearestNeighbor.h"
#include "ExtraMileage.h"

#include <pthread.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


/* 
 * Instead of using time() try with:
 *
    struct timespec start, finish;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &finish);
    double elapsed = ((finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec) / 1000000000.0); // nsed are useless if we are counting time in minutes
    return elapsed; // THIS IS IN SECONDS
 * 
 * This should be more robust(according to a guy on stackoverflow)
*/

/* Checks wheter the time passed since the initialization of start has passed timeLimit.
    Returns 1 in that case, 0 otherwise.*/
/*static inline int checkTime(struct timespec start, double timeLimit);

Solution VariableNeighborhood(Instance *inst, int configuration)
{
    // time limit management
    double placeholderTime = 100.0;
    struct timespec start;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &start);

    Solution sol;

    // check if a seed for random has been passed as argument
    if(inst->params.randomSeed != -1) srand(inst->params.randomSeed);
    else throwError(inst, &sol, "VariableNeighborhood: random seed has not been passed as argument");
    
    if(configuration != 0 && configuration != 1) throwError(inst, &sol, "VariableNeighborhood: incorrect argument for configuration");
    else if(configuration == 0) // Find local minimum with Nearest Neighbour
    {
        // Compute a solution with Nearest Neighbour and optimize it with 2-opt
        sol = NearestNeighbour(inst);
        _2optBestFix(inst);
        while(checkTime(start, placeholderTime) == 0)
        {
            changeNeighborhood(inst);
            _2optBestFix(inst);
            
        }

    }else // Find local minimum with Extra Milage
    {
        // Compute a solution with Extra Mileage and optimize it with 2-opt
        sol = solveExtraMileage(inst);
        _2optBestFix(inst);
        while(checkTime(start, placeholderTime) == 0)
        {
            _2optBestFix(inst);
            
        }
    }
    return sol;
}

static inline void changeNeighborhood(Solution *sol)
{

}

static inline int checkTime(struct timespec start, double timeLimit)
{
    struct timespec currentTime;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currentTime);
    double elapsed = ((currentTime.tv_sec - start.tv_sec) + (currentTime.tv_nsec - start.tv_nsec) / 1000000000.0); // nsed are useless if we are counting time in minutes
    
    if(elapsed < timeLimit) return 0;
    else return 1;
}
*/