#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance d; initInstance(&d);
    parseArgs(&d, argc, argv);
    readFile(&d);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", d.params.inputFile);

    double computeMatrixTime = computeDistanceMatrix(&d);
    
    LOG(LOG_LVL_LOG, "Distance Matrix Done in %.3e seconds", computeMatrixTime);

    // WORKS WITH THIS COMMENTED
    NearestNeighbour(&d);
    _2optBestFix(&d);
    
    for(int i = 0; i < d.nodesCount; i++) LOG(LOG_LVL_EVERYTHING, "Node %d in solution: %d", i, d.solution.bestSolution[i]);

    //saveSolution(&d);
    plotSolution(&d, "640,480", "green", "black", 1);    
    //###########################################

    freeInstance(&d);

    return EXIT_SUCCESS;
}