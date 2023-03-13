#include "tsp.h"

//int distMatIntegerW(Instance *d);

// method to compute edges weights with multiple threads
static void * fillMatrix(void * d);

int computeDistanceMatrix(Instance *d)
{
    /*  uncomment method on top also 
    if (Type requires so)
       return distMatIntegerW;
    */

    // allocate memory
    d->edgeCost = malloc(d->nodesCount * d->nodesCount * sizeof(double) + 4 * sizeof(double)); // using aligned alloc to use load() instead of uload() will make things messier

    // initialize mutex
    pthread_mutex_init(&d->global.mutex, NULL);

}

static void * fillMatrix(void* arg)
{

}