#include "TspCplex.h"
#include "Tsp.h"

#include <cplex.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time



void HardFixing(Solution *sol, double fixingAmount, enum HardFixPolicy policy, double tlim)
{
    struct timespec currT;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    double startTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
	double currentTimeSec = startTimeSec;

    Instance *inst = sol->instance;
    size_t n = inst->nNodes;

    if (inst->params.logLevel >= LOG_LVL_DEBUG)
        checkSolution(sol);

    CplexData cpx = initCplexData(inst);
    int ncols = CPXgetnumcols(cpx.env, cpx.lp);

    // allocate necessary memory
    double *xstar = malloc(ncols * sizeof(double));
    int *indexes = malloc(ncols * sizeof(int));
    SubtoursData sub = {
        .subtoursCount = 0,
        .successors = malloc(n * sizeof(int)),
        .subtoursMap = malloc(n * sizeof(int))
    };
    int *bestSuccessorsSol = malloc(n * sizeof(int));
    double bestSuccCost = INFINITY;

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
    currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;

    while (currentTimeSec > startTimeSec + tlim)
    {
        switch (policy)
        {
        case HARDFIX_POLICY_RANDOM:
            break;
        case HARDFIX_POLICY_SMALLEST:
            break;
        case HARDFIX_POLICY_PROBABILITY:
            break;
        case HARDFIX_POLICY_MIXED:
            break;
        }

        double remainingTime = startTimeSec + tlim - currentTimeSec;
        CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, remainingTime);

        clock_gettime(_POSIX_MONOTONIC_CLOCK, &currT);
        currentTimeSec = (double)currT.tv_sec + (double)currT.tv_nsec / 1000000000.0;
    }
    destroyCplexData(&cpx);
    free(xstar); free(indexes); free(sub.successors); free(sub.subtoursMap);
    
    if (bestSuccCost < sol->cost)
        cvtSuccessorsToSolution(bestSuccessorsSol, sol);
    else
        LOG(LOG_LVL_WARNING, "HardFixing: Solution could not be optimized any further");

    free(bestSuccessorsSol);
}