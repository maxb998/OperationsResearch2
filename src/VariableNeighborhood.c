#include "tsp.h"

double VariableNeighborhood(Instance *inst, int configuration)
{
    intmax_t placeholder = 1000000;
    // check if a seed for random has been passed as argument
    if(inst->params.randomSeed != -1) srand(inst->params.randomSeed);
    else throwError(inst, "VariableNeighborhood: random seed has not been passed as argument");
    
    if(configuration != 0 && configuration != 1) throwError(inst, "VariableNeighborhood: incorrect argument for configuration");
    else if(configuration == 0) // Find local minimum with Nearest Neighbour
    {
        // Initialization of the time for the time limit
        time_t startingTime = time(NULL);

        // Compute a solution with Nearest Neighbour and optimize it with 2-opt
        NearestNeighbour(inst);
        _2optBestFix(inst);
        while((intmax_t)time(NULL) - startingTime < placeholder)
        {
            changeNeighborhood(inst);
            _2optBestFix(inst);
            
        }

    }else                       // Find local minimum with Extra Milage
    {
        // Initialization of the time for the time limit
        time_t startingTime = time(NULL);

        // Compute a solution with Extra Mileage and optimize it with 2-opt
        solveExtraMileage(inst);
        _2optBestFix(inst);
        while((intmax_t)time(NULL) - startingTime < placeholder)
        {
            _2optBestFix(inst);
            
        }
    }

    return 0;
}

static inline void changeNeighborhood(Instance *inst)
{

}
