#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance d; initInstance(&d);
    parseArgs(&d, argc, argv);
    readFile(&d);

    //LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", d.params.inputFile);

    //double computeMatrixTime = computeDistanceMatrix(&d);
    
    //LOG(LOG_LVL_LOG, "Distance Matrix Done in %.3e seconds", computeMatrixTime);

    //printDistanceMatrix(&d, 0);

    //###########################################
    // Use best solution to see if it plots correctly
    d.solution.bestSolution = malloc(d.nodesCount * sizeof(int));
    
    // real ugly
    int optSolution[48] = { 0,7,37,30,43,17,6,27,5,36,18,26,16,42,29,35,45,32,19,46,20,31,38,47,4,41,23,9,44,34,3,25,1,28,33,40,15,21,2,22,13,24,12,10,11,14,39,8 };

    memcpy(d.solution.bestSolution, optSolution, d.nodesCount * sizeof(int));

    //saveSolution(&d);
    plotSolution(&d, "1920,1080", "green", "black", 2);    
    //###########################################

    freeInstance(&d);

    return EXIT_SUCCESS;
}