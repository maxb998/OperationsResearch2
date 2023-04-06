#include "TspUtilities.h"
#include "TspFileUtils.h"

#include "CostMatrix.h"

#include "NearestNeighbor.h"
#include "ExtraMileage.h"
#include "2Opt.h"

int main (int argc, char *argv[])
{
    Instance inst = newInstance();
    parseArgs(&inst, argc, argv);
    readFile(&inst);

    LOG (LOG_LVL_LOG, "file %s has been loaded succesfully", inst.params.inputFile);

    double computeMatrixTime = computeCostMatrix(&inst);


    
    
    LOG(LOG_LVL_LOG, "Distance Matrix done in %lf seconds", computeMatrixTime);
    
    
    Solution nn = NearestNeighbor(&inst);
    LOG(LOG_LVL_LOG, "Nearest Neighbour finished in %lf seconds", nn.execTime);

    double _2optTime = apply2OptBestFix(&nn, _2OPT_AVX_ST);
    LOG(LOG_LVL_LOG, "2-Opt finished in %lf seconds", _2optTime);

    Solution em = ExtraMileage(&inst, EM_INIT_EXTREMES);
    
    plotSolution(&nn, "1366,768", "green", "black", 1);
    plotSolution(&em, "1920,1080", "green", "black", 1);

    destroySolution(&em);
    destroyInstance(&inst);


    return EXIT_SUCCESS;
}