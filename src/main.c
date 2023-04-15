#include "TspUtilities.h"
#include "TspIOUtils.h"

#include "CostMatrix.h"

#include "ArgParser.h"
#include "NearestNeighbor.h"
#include "ExtraMileage.h"
#include "2Opt.h"

int main (int argc, char *argv[])
{
    Instance inst = newInstance();
    argParse(&inst, argc, argv);
    readFile(&inst);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", inst.params.inputFile);

    double computeMatrixTime = computeCostMatrix(&inst);


    
    
    LOG(LOG_LVL_LOG, "Distance Matrix done in %lf seconds", computeMatrixTime);
    
    
    Solution nn = NearestNeighbor(&inst);
    LOG(LOG_LVL_LOG, "Nearest Neighbor finished in %lf seconds. Solution cost is %lf", nn.execTime, nn.bestCost);

    double _2optTime = apply2OptBestFix(&nn, _2OPT_AVX_ST);
    LOG(LOG_LVL_LOG, "2-Opt finished optimizing Nearest Neighbor in %lf seconds. Solution cost is %lf", _2optTime, nn.bestCost);

    /*Solution em = ExtraMileage(&inst, EM_OPTION_AVX, EM_INIT_RANDOM);
    LOG(LOG_LVL_LOG, "Extra Mileage finished in %lf seconds. Solution cost is %lf", em.execTime, em.bestCost);
    LOG(LOG_LVL_LOG, "Cost of Extra Mileage Solution is %lf", computeSolutionCostVectorizedFloat(&em));*/
    
    plotSolution(&nn, "3600,2000", "green", "black", 1, 0);
    //plotSolution(&em, "1920,1080", "green", "black", 1, 0);

    destroySolution(&nn);
    //destroySolution(&em);
    destroyInstance(&inst);


    return EXIT_SUCCESS;
}