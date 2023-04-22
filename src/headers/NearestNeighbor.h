#ifndef TSP_HEURISTICS_H
#define TSP_HEURISTICS_H

#include "TspBase.h"
#include <limits.h>

Solution NearestNeighbor(Instance *inst, enum NNFirstNodeOptions startOption, double timeLimit, int useThreads);

#endif // TSP_HEURISTICS_H