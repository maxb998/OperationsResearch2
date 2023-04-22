#ifndef EXTRA_MILEAGE
#define EXTRA_MILEAGE

#include "TspBase.h"

enum EMOptions {
    EM_OPTION_AVX,
    EM_OPTION_BASE,
    EM_OPTION_USE_COST_MATRIX
};

Solution ExtraMileage(Instance *inst, enum EMOptions emOpt, enum EMInitType emInitType);

/* 
 * Function that applies Extra Mileage heuristic on a solution that is already not empty.
 * ARGUMENTS:
 *  sol -> solution element to which apply extra mileage
 *  nCovered -> number of elements representing the first part of the tour inside sol. The function will add all remaning nodes
 *  emOpt -> specify the way thje algorithm should operate (zero is the fastest)
 * RETURN VALUE: Number in seconds representing execution time
*/
double applyExtraMileage(Solution *sol, size_t nCovered, enum EMOptions emOpt);

#endif // EXTRA_MILEAGE