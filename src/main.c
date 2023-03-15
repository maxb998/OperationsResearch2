#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance d; initInstance(&d);
    parseArgs(&d, argc, argv);
    readFile(&d);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", d.params.inputFile);

    computeSquaredDistanceMatrix(&d);
    
    LOG(LOG_LVL_LOG, "Distance Matrix Done");

    //printDistanceMatrix(&d, 1);

    //###########################################
    // CREATE A FAKE SOLUTION USING THE POINTS IN THEIR DEFAULT ORDER
    int fakeSolution[48];
    for (size_t i = 0; i < 48; i++)
    {
        fakeSolution[i] = i;
    }
    d.solution.bestSolution = fakeSolution;
    //saveSolution(&d);
    //plotSolution(&d);    
    //###########################################

    freeInstance(&d);

    return EXIT_SUCCESS;
}