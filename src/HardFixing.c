#include "TspCplex.h"

#include <cplex.h>

Solution HardFixing(Solution *sol, double tlim)
{
    Instance *inst = sol->instance;

    CplexData cpx = initCplexData(inst);




    destroyCplexData(&cpx);
}