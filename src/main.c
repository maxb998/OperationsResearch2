#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance inst = newInstance();
    parseArgs(&inst, argc, argv);
    readFile(&inst);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", inst.params.inputFile);

    //double computeMatrixTime = computeDistanceMatrix(&d);
    
    //LOG(LOG_LVL_LOG, "Distance Matrix Done in %.3e seconds", computeMatrixTime);

    // WORKS WITH THIS COMMENTED
    Solution nn = NearestNeighbour(&inst);
    //_2optBestFix(&d);
    
    //saveSolution(&nn);

    LOG(LOG_LVL_WARNING, "Total cost of Solution recomputed is %f", computeSquaredCost_VEC(&nn));

    //for(int i = 0; i < d.nNodes; i++) LOG(LOG_LVL_EVERYTHING, "Node %d in solution: %d", i, d.solution.bestSolution[i]);

    //saveSolution(&d);
    plotSolution(&nn, "1920,1080", "green", "black", 1);    
    //###########################################
    destroySolution(&nn);
    destroyInstance(&inst);


    return EXIT_SUCCESS;
}