#include "tsp.h"

double VariableNeighborhood(Instance *d, int configuration)
{
    intmax_t placeholder = 1000000;
    // check if a seed for random has been passed as argument
    if(d->params.randomSeed != -1) srand(d->params.randomSeed);
    else throwError(d, "VariableNeighborhood: random seed has not been passed as argument");
    
    if(configuration != 0 && configuration != 1) throwError(d, "VariableNeighborhood: incorrect argument for configuration");
    else if(configuration == 0) // Find local minimum with Nearest Neighbour
    {
        // Initialization of the time for the time limit
        time_t startingTime = time(NULL);

        
        NearestNeighbour(d);

    }else                       // Find local minimum with Extra Milage
    {

    }

    return 0;
}