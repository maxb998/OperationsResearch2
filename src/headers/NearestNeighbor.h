#ifndef TSP_HEURISTICS_H
#define TSP_HEURISTICS_H

#include "TspBase.h"
#include <limits.h>

#define NN_GRASP_COEFF (int)((double)RAND_MAX/0.8)

Solution NearestNeighbour(Instance *inst);

#endif // TSP_HEURISTICS_H