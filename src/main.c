#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance d; d.params.useAVX = 0;    // avx disabled by default
    parseArgs(&d, argc, argv);
    readFile(&d);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", d.params.inputFile);

    //###########################################
    // CREATE A FAKE SOLUTION USING THE POINTS IN THEIR DEFAULT ORDER
    int fakeSolution[48];
    for (size_t i = 0; i < 48; i++)
    {
        fakeSolution[i] = i;
    }
    d.solution.bestSolution = fakeSolution;
    saveSolution(&d);
    plotSolution(&d);    
    //###########################################

    return EXIT_SUCCESS;
}