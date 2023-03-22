#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance d; initInstance(&d);
    parseArgs(&d, argc, argv);
    readFile(&d);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", d.params.inputFile);

    double computeMatrixTime = computeDistanceMatrix(&d);
    
    LOG(LOG_LVL_LOG, "Distance Matrix Done in %.3e seconds", computeMatrixTime);

    printDistanceMatrix(&d, 0);

    //###########################################
    // CREATE A FAKE SOLUTION USING THE POINTS IN THEIR DEFAULT ORDER
    int *fakeSolution = malloc(48 * sizeof(int));
    for (int i = 0; i < 48; i++)
    {
        fakeSolution[i] = i+1;
    }
    fakeSolution[47] = 0;
    d.solution.bestSolution = fakeSolution;
    //saveSolution(&d);
    //plotSolution(&d);    
    //###########################################

    freeInstance(&d);

    return EXIT_SUCCESS;
}