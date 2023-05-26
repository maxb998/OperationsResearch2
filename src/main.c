#include "Tsp.h"
#include "TspCplex.h"

#include <stdio.h>
#include <time.h>


static Solution runNearestNeighbor(Instance *inst);
static Solution runExtraMileage(Instance *inst);
static Solution runVariableNeighborhoodSearch(Instance *inst);

static Solution runBenders(Instance *inst);
static Solution runBranchAndCut(Instance *inst);
static Solution runHardFixing(Instance *inst);


static void run2Opt(Solution *sol);


int main (int argc, char *argv[])
{
    Instance inst = newInstance();
    argParse(&inst, argc, argv);

    if (inst.params.randomSeed != -1)
        srand(inst.params.randomSeed);
    else
        srand(time(NULL));

    printInfo(&inst);

    double fileReadTime = readFile(&inst);
    LOG (LOG_LVL_NOTICE, "file %s has been loaded succesfully in %lf milliseconds", inst.params.inputFile, fileReadTime * 1000.);
    
    double computeMatrixTime = computeCostMatrix(&inst);
    LOG(LOG_LVL_NOTICE, "Distance Matrix done in %lf seconds", computeMatrixTime);
    
    // initializing pointers to null to avoid possible errors on destruction of sol at the end of main
    Solution sol = { .indexPath = NULL, .X = NULL, .Y = NULL };

    switch (inst.params.mode)
    {
    case MODE_NONE:
        // should never enter here
        throwError(&inst, NULL, "mode is not set even after checking in argParse");
        break;

    case MODE_NN:
        sol = runNearestNeighbor(&inst);
        if (inst.params.use2OptFlag)
            run2Opt(&sol);
        break;

    case MODE_EM:
        sol = runExtraMileage(&inst);
        if (inst.params.use2OptFlag)
            run2Opt(&sol);
        break;

    case MODE_VNS:
        sol = runVariableNeighborhoodSearch(&inst);
        break;

    case MODE_BENDERS:
        sol = runBenders(&inst);
        break;

    case MODE_BRANCH_CUT:
        sol = runBranchAndCut(&inst);
        break;

    case MODE_HARDFIX:
        sol = runHardFixing(&inst);
        break;
    }


    if (inst.params.showPlotFlag)
        plotSolution(&sol, "800,600", "green", "black", 1, 0);
    
    if (inst.params.saveFlag)
        saveSolution(&sol, argc, argv);

    destroySolution(&sol);
    destroyInstance(&inst);

    printf("\n");

    return EXIT_SUCCESS;
}


static Solution runNearestNeighbor(Instance *inst)
{
    printf("\nNearest Neighbor starting...\n");

    Solution nn = NearestNeighbor(inst, inst->params.nnFirstNodeOption, inst->params.tlim, 1);

    printf("Nearest Neighbot finished in %lf second\n", nn.execTime);
    printf("Cost = %lf\n", nn.cost);

    return nn;
}

static Solution runExtraMileage(Instance *inst)
{
    printf("\nExtra Mileage starting...\n");

    Solution em = ExtraMileage(inst, 0, inst->params.emInitOption);

    printf("Extra Mileage finished in %lf seconds\n", em.execTime);
    printf("Cost = %lf\n", em.cost);

    return em;
}

static Solution runVariableNeighborhoodSearch(Instance *inst)
{
    Solution sol = VariableNeighborhood(inst, inst->params.vnsInitOption);
    return sol;
}

static Solution runBenders(Instance *inst)
{
    printf("##############################################################################################################################\n");
    printf("Benders starting...\n");
    printf("##############################################################################################################################\n");

    Solution sol = benders(inst, inst->params.tlim);

    printf("##############################################################################################################################\n");
    printf("Benders finished in %lf seconds\n", sol.execTime);
    printf("Cost = %lf\n", sol.cost);
    printf("##############################################################################################################################\n");

    return sol;
}

static Solution runBranchAndCut(Instance *inst)
{
    printf("Running TEMPORARELY extra mileage to get a warm start solution\n");
    Solution warmStart = ExtraMileage(inst, 0, EM_INIT_EXTREMES);
    apply2OptBestFix(&warmStart, 0);

    printf("##############################################################################################################################\n");
    printf("Branch & Cut starting...\n");
    printf("##############################################################################################################################\n");

    Solution sol = BranchAndCut(inst, inst->params.tlim, &warmStart);

    destroySolution(&warmStart);

    printf("##############################################################################################################################\n");
    printf("Branch & Cut finished in %lf seconds\n", sol.execTime);
    printf("Cost = %lf\n", sol.cost);
    printf("##############################################################################################################################\n");

    return sol;
}

static Solution runHardFixing(Instance *inst)
{
    printf("Running TEMPORARELY extra mileage to get a warm start solution\n");
    Solution sol = ExtraMileage(inst, 0, EM_INIT_EXTREMES);
    apply2OptBestFix(&sol, 0);

    printf("##############################################################################################################################\n");
    printf("Hard Fixing Starting...\n");
    printf("##############################################################################################################################\n");

    HardFixing(&sol, 0.5, HARDFIX_POLICY_RANDOM, inst->params.tlim);

    printf("##############################################################################################################################\n");
    printf("Hard Fixing finished in %lf seconds\n", sol.execTime);
    printf("Cost = %lf\n", sol.cost);
    printf("##############################################################################################################################\n");

    return sol;
}

static void run2Opt(Solution *sol)
{
    printf("\n2Opt starting...\n");

    double optTime = apply2OptBestFix(sol, _2OPT_AVX_ST);
    sol->execTime += optTime;

    printf("2Opt finished in %lf seconds\n", optTime);
    printf("Cost = %lf\n", sol->cost);
}