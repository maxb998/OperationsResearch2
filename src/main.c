#include "Tsp.h"
#include "TspCplex.h"

#include <stdio.h>
#include <time.h>

#define SEPARATOR_STR "##############################################################################################################################\n"

#define METAHEUR_INIT_RATIO 0.05
#define MATHEUR_METAHEUR_INIT_RATIO 0.1
#define MATHEUR_INIT_RATIO 0.01

static Solution runHeuristic(Instance *inst, enum Mode mode, double tlim);
static void runMetaheuristic(Solution *sol, enum Mode mode, double tlim);
static void runExactSolver(Solution *sol, enum Mode mode, double tlim);
static void runMatheuristic (Solution *sol, enum Mode mode, double tlim);

static void run2Opt(Solution *sol);
static void run3Opt(Solution *sol);


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

    if (inst.params.edgeWeightType > 4)
        throwError("Edge Weight Type is not supported. This solver only support: EUC_2D, MAN_2D, MAX_2D, CEIL_2D, ATT");
    
    #if (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
        printf(SEPARATOR_STR);
        double computeMatrixTime = computeCostMatrix(&inst);
        LOG(LOG_LVL_NOTICE, "Distance Matrix done in %lf seconds", computeMatrixTime);
        printf(SEPARATOR_STR"\n");
    #elif ((COMPUTATION_TYPE == COMPUTE_OPTION_USE_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_USE_BASE))
    #endif

    // initializing pointers to null to avoid possible errors on destruction of sol at the end of main
    Solution sol = { .instance=&inst, .indexPath = NULL };

    double tlim = inst.params.tlim;
    enum Mode m = inst.params.mode;

    if (m & (MODE_NN | MODE_EM))
    {
        sol = runHeuristic(&inst, m, tlim);
    }
    else if (m & (MODE_TABU | MODE_VNS | MODE_ANNEALING))
    {
        sol = runHeuristic(&inst, inst.params.metaheurInitMode, tlim * METAHEUR_INIT_RATIO);
        run2Opt(&sol);
        runMetaheuristic(&sol, m, tlim - sol.execTime);
    }
    else if (m & MODE_GENETIC)
    {
        runMetaheuristic(&sol, MODE_GENETIC, tlim);
    }
    else
    {
        if ((inst.params.cplexWarmStart || (m & (MODE_HARDFIX | MODE_LOCAL_BRANCHING))) && (inst.params.matheurInitMode & (MODE_NN | MODE_EM)))
        {
            sol = runHeuristic(&inst, inst.params.matheurInitMode, tlim * MATHEUR_INIT_RATIO);
            run2Opt(&sol);
        }
        else if ((inst.params.cplexWarmStart || (m & (MODE_HARDFIX | MODE_LOCAL_BRANCHING))) && (inst.params.matheurInitMode & (MODE_TABU | MODE_VNS | MODE_ANNEALING)))
        {
            sol = runHeuristic(&inst, inst.params.metaheurInitMode, (tlim * MATHEUR_INIT_RATIO) * MATHEUR_METAHEUR_INIT_RATIO);
            run2Opt(&sol);
            runMetaheuristic(&sol, inst.params.matheurInitMode, (tlim * MATHEUR_INIT_RATIO) * (1 - MATHEUR_METAHEUR_INIT_RATIO));
        }
        else if (inst.params.cplexWarmStart && (inst.params.matheurInitMode & MODE_GENETIC))
            runMetaheuristic(&sol, MODE_GENETIC, tlim * MATHEUR_INIT_RATIO);
        else if (!inst.params.cplexWarmStart)
            sol = newSolution(&inst);

        if (m & (MODE_BENDERS | MODE_BRANCH_CUT))
            runExactSolver(&sol, m, tlim - sol.execTime);
        else // if (m & (MODE_HARDFIX | MODE_LOCAL_BRANCHING))
            runMatheuristic(&sol, m, tlim - sol.execTime);
    }

    if (inst.params.use2Opt)
        run2Opt(&sol);
    
    if (inst.params.use3Opt)
        run3Opt(&sol);

    if (inst.params.showPlot)
        plotSolution(&sol, "1920,1080", "green", "black", 1, false);
    
    if (inst.params.saveSolution)
        saveSolution(&sol, argc, argv);
    
    printf(SEPARATOR_STR);
    printf("Total runtime = %lf seconds\n", sol.execTime);
    printf("Final cost = %lf\n", cvtCost2Double(sol.cost));
    printf(SEPARATOR_STR);

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

        sol = NearestNeighbor(inst, tlim);

        printf("Nearest Neighbor finished in %lf second\n", sol.execTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol.cost));
        break;

    case MODE_EM:
        printf("Extra Mileage starting...\n");

        sol = ExtraMileage(inst, tlim);

        printf("Extra Mileage finished in %lf seconds\n", sol.execTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol.cost));
        break;

    default: throwError("runHeuristic: specified mode must be in [MODE_NN, MODE_EM]"); break;
    }

    printf(SEPARATOR_STR"\n");

    return sol;
}

static void runMetaheuristic(Solution *sol, enum Mode mode, double tlim)
{
    double startTime = sol->execTime;

    printf(SEPARATOR_STR);

    switch (mode)
    {
    case MODE_TABU:
        printf("Tabu Search Starting...\n");

        TabuSearch(sol, tlim);

        printf("Tabu Search finished in %lf second\n", sol->execTime - startTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;
    case MODE_VNS:
        printf("Variable Neighborhood Search Starting...\n");

        VariableNeighborhoodSearch(sol, tlim);

        printf("Variable Neighborhood Search finished in %lf second\n", sol->execTime - startTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;
    case MODE_ANNEALING:
        printf("Simulated Annealing Starting...\n");

        SimulatedAnnealing(sol, tlim);

        printf("Simulated Annealing finished in %lf second\n", sol->execTime - startTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;
    case MODE_GENETIC:
        printf("Genetic Search Starting...\n");

        *sol = GeneticAlgorithm(sol->instance, tlim);

        printf("Genetic Search finished in %lf second\n", sol->execTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;


    default: throwError("runMetaheuristic: specified mode must be in [MODE_TABU, MODE_VNS, MODE_ANNEALING, MODE_GENETIC]"); break;
    }

    printf(SEPARATOR_STR"\n");
}

static void runExactSolver(Solution *sol, enum Mode mode, double tlim)
{
    double startTime = sol->execTime;

    printf(SEPARATOR_STR);

    switch (mode)
    {
    case MODE_BENDERS:
        printf("Benders Starting...\n");

        benders(sol, tlim);

        printf("Benders finished in %lf second\n", sol->execTime - startTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;
    case MODE_BRANCH_CUT:
        printf("Branch & Cut Starting...\n");

        BranchAndCut(sol, tlim);

        printf("Branch & Cut finished in %lf second\n", sol->execTime - startTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;

    default: throwError("runMetaheuristic: specified mode must be in [MODE_BENDERS, MODE_BRANCH_CUT]"); break;
    }

    printf(SEPARATOR_STR"\n");
}

static void runMatheuristic (Solution *sol, enum Mode mode, double tlim)
{
    double startTime = sol->execTime;

    printf(SEPARATOR_STR);

    switch (mode)
    {
    case MODE_HARDFIX:
        printf("Hard Fixing Starting...\n");

        HardFixing(sol, tlim);

        printf("Hard Fixing finished in %lf second\n", sol->execTime - startTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;
    case MODE_LOCAL_BRANCHING:
        printf("Local Branching Starting...\n");

        LocalBranching(sol, tlim);

        printf("Local Branching finished in %lf second\n", sol->execTime - startTime);
        printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
        break;

    default: throwError("runMetaheuristic: specified mode must be in [MODE_HARDFIX, MODE_LOCAL_BRANCHING]"); break;
    }

    printf(SEPARATOR_STR"\n");
}

static void run2Opt(Solution *sol)
{
    printf(SEPARATOR_STR);
    printf("2Opt starting...\n");

    double startTime = sol->execTime;
    set2OptPerformanceBenchmarkLogMT(true);
    apply2OptBestFixMT(sol);
    set2OptPerformanceBenchmarkLogMT(false);

    printf("2Opt finished in %lf seconds\n", sol->execTime - startTime);
    printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
    printf(SEPARATOR_STR"\n");
}

static void run3Opt(Solution *sol)
{
    printf(SEPARATOR_STR);
    printf("3Opt starting...\n");

    double startTime = sol->execTime;
    set3OptPerformanceBenchmarkLogMT(true);
    apply3OptBestFixMT(sol);
    set3OptPerformanceBenchmarkLogMT(false);

    printf("3Opt finished in %lf seconds\n", sol->execTime - startTime);
    printf("Solution Cost = %lf\n", cvtCost2Double(sol->cost));
    printf(SEPARATOR_STR"\n");
}

