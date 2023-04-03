#include "TspUtilities.h"
#include "TspFileUtils.h"

#include "CostMatrix.h"

#include "NearestNeighbor.h"
#include "2Opt.h"

int main (int argc, char *argv[])
{
    Instance inst = newInstance();
    parseArgs(&inst, argc, argv);
    readFile(&inst);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", inst.params.inputFile);

    double computeMatrixTime = computeCostMatrix(&inst);
    
    LOG(LOG_LVL_LOG, "Distance Matrix done in %lf seconds", computeMatrixTime);
    
    //printCostMatrix(&inst);

    // WORKS WITH THIS COMMENTED
    Solution nn = NearestNeighbour(&inst);
    LOG(LOG_LVL_LOG, "Nearest Neighbour finished in %lf seconds", nn.execTime);
    double _2optTime = apply2OptBestFix(&nn, _2OPT_BASE);
    LOG(LOG_LVL_LOG, "2-Opt finished in %lf seconds", _2optTime);
    
    //saveSolution(&nn);

    //LOG(LOG_LVL_WARNING, "Total cost of Solution recomputed is %f", computeSolutionCostVectorizedDouble(&nn));
    //LOG(LOG_LVL_WARNING, "Total cost of Solution recomputed is %lf", computeSolutionCost(&nn));

    //for(int i = 0; i < d.nNodes; i++) LOG(LOG_LVL_EVERYTHING, "Node %d in solution: %d", i, d.solution.bestSolution[i]);

    //saveSolution(&d);
    plotSolution(&nn, "1920,1080", "green", "black", 1);    
    //###########################################
    //destroySolution(&nn);
    destroyInstance(&inst);


    return EXIT_SUCCESS;
}