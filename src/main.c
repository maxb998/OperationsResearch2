#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance d; d.params.useAVX = 0;    // avx disabled by default
    parseArgs(&d, argc, argv);
    readFile(&d);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", d.params.inputFile);

    computeSquaredDistanceMatrix(&d);
    
    LOG(LOG_LVL_LOG, "Distance Matrix Done");

    printDistanceMatrix(&d, 1);

    return EXIT_SUCCESS;
}