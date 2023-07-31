#include "Tsp.h"
#include "TspCplex.h"

#include <stdio.h>
#include <time.h>

#define SEPARATOR_STR "##############################################################################################################################\n"

#define METAHEUR_INIT_RATIO 1/10
#define MATHEUR_INIT_RATIO 1/10

#define MATH

static Solution runHeuristic(Instance *inst, enum Mode mode, double tlim);
static void runMetaheuristic(Solution *sol, enum Mode mode, double tlim);
static void runExactSolver(Solution *sol, enum Mode mode, double tlim);
static void runMatheuristic (Solution *sol, enum Mode mode, double tlim);

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
    printf("\n");
    
    if (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
    {
        printf(SEPARATOR_STR);
        double computeMatrixTime = computeCostMatrix(&inst);
        LOG(LOG_LVL_NOTICE, "Distance Matrix done in %lf seconds", computeMatrixTime);
        printf(SEPARATOR_STR"\n");
    }

    // initializing pointers to null to avoid possible errors on destruction of sol at the end of main
    Solution sol = { .indexPath = NULL };

    double tlim = inst.params.tlim;
    enum Mode m = inst.params.mode;

    if ((m == MODE_NN) || (m == MODE_EM))
    {
        sol = runHeuristic(&inst, m, tlim);
    }
    else if ((m >= MODE_TABU) && (m <= MODE_GENETIC))
    {
        sol = runHeuristic(&inst, inst.params.metaheurInitMode, tlim * METAHEUR_INIT_RATIO );
        run2Opt(&sol);
        runMetaheuristic(&sol, m, tlim * (1 - METAHEUR_INIT_RATIO));
    }
    else
    {
        enum Mode init;
        if (m <= MODE_BRANCH_CUT)
            init = inst.params.warmStartMode;
        else
            init = inst.params.matheurInitMode;

        double remainingTime = tlim;
        if (init <= MODE_EM)
        {
            sol = runHeuristic(&inst, init, tlim * MATHEUR_INIT_RATIO);
            run2Opt(&sol);
        }
        else 
        {
            sol = runHeuristic(&inst, inst.params.metaheurInitMode, tlim * METAHEUR_INIT_RATIO);
            runMetaheuristic(&sol, m, tlim * (1 - METAHEUR_INIT_RATIO));
            run2Opt(&sol);
            remainingTime -= tlim * METAHEUR_INIT_RATIO;
        }
        remainingTime -= tlim * MATHEUR_INIT_RATIO;

        if (m <= MODE_BRANCH_CUT)
            runExactSolver(&sol, m, remainingTime);
        else
            runMatheuristic(&sol, m, remainingTime);
    }

    if (inst.params.use2Opt)
        run2Opt(&sol);

    if (inst.params.showPlot)
        plotSolution(&sol, "1920,1080", "green", "black", 1, false);
    
    if (inst.params.saveSolution)
        saveSolution(&sol, argc, argv);

    destroySolution(&sol);
    destroyInstance(&inst);

    printf("\n");

    return EXIT_SUCCESS;
}

static Solution runHeuristic(Instance *inst, enum Mode mode, double tlim)
{
    Solution sol = {0};

    printf(SEPARATOR_STR);

    switch (mode)
    {
    case MODE_NN:
        printf("Nearest Neighbor Starting...\n");

        sol = NearestNeighbor(inst, inst->params.nnFirstNodeOption, tlim, 0);

        printf("Nearest Neighbor finished in %lf second\n", sol.execTime);
        printf("Solution Cost = %lf\n", sol.cost);
        break;

    case MODE_EM:
        printf("Extra Mileage starting...\n");

        sol = ExtraMileage(inst, inst->params.emInitOption, tlim, 0);

        printf("Extra Mileage finished in %lf seconds\n", sol.execTime);
        printf("Solution Cost = %lf\n", sol.cost);
        break;

    default: throwError(inst, NULL, "runHeuristic: specified mode must be in [MODE_NN, MODE_EM]"); break;
    }

    printf(SEPARATOR_STR"\n");

    return sol;
}

static void runMetaheuristic(Solution *sol, enum Mode mode, double tlim)
{
    printf(SEPARATOR_STR);

    switch (mode)
    {
    case MODE_TABU:
        printf("Tabu Search Starting...\n");

        printf("Tabu Search finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;
    case MODE_VNS:
        printf("Variable Neighborhood Search Starting...\n");

        VariableNeighborhoodSearch(sol, tlim, 0);

        printf("Variable Neighborhood Search finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;
    case MODE_ANNEALING:
        printf("Simulated Annealing Starting...\n");

        printf("Simulated Annealing finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;
    case MODE_GENETIC:
        printf("Genetic Search Starting...\n");

        printf("Genetic Search finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;

    default: throwError(sol->instance, sol, "runMetaheuristic: specified mode must be in [MODE_TABU, MODE_VNS, MODE_ANNEALING, MODE_GENETIC]"); break;
    }

    printf(SEPARATOR_STR"\n");
}

static void runExactSolver(Solution *sol, enum Mode mode, double tlim)
{
    printf(SEPARATOR_STR);

    switch (mode)
    {
    case MODE_BENDERS:
        printf("Benders Starting...\n");

        benders(sol, tlim);

        printf("Benders finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;
    case MODE_BRANCH_CUT:
        printf("Branch & Cut Starting...\n");

        BranchAndCut(sol, tlim);

        printf("Branch & Cut finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;

    default: throwError(sol->instance, sol, "runMetaheuristic: specified mode must be in [MODE_BENDERS, MODE_BRANCH_CUT]"); break;
    }

    printf(SEPARATOR_STR"\n");
}

static void runMatheuristic (Solution *sol, enum Mode mode, double tlim)
{
    printf(SEPARATOR_STR);

    switch (mode)
    {
    case MODE_HARDFIX:
        printf("Hard Fixing Starting...\n");

        HardFixing(sol, 0.5, sol->instance->params.hardFixPolicy, tlim);

        printf("Hard Fixing finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;
    case MODE_LOCAL_BRANCHING:
        printf("Local Branching Starting...\n");

        printf("Local Branching finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", sol->cost);
        break;

    default: throwError(sol->instance, sol, "runMetaheuristic: specified mode must be in [MODE_HARDFIX, MODE_LOCAL_BRANCHING]"); break;
    }

    printf(SEPARATOR_STR"\n");
}

static void run2Opt(Solution *sol)
{
    printf(SEPARATOR_STR);
    printf("2Opt starting...\n");

    set2OptPerformanceBenchmarkLog(true);
    apply2OptBestFix(sol);
    set2OptPerformanceBenchmarkLog(false);

    printf("2Opt finished in %lf seconds\n", sol->execTime);
    printf("Solution Cost = %lf\n", sol->cost);
    printf(SEPARATOR_STR"\n");
}

