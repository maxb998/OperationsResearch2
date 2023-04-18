#ifndef TSP_HEURISTICS_H
#define TSP_HEURISTICS_H

#include "TspBase.h"
#include <limits.h>

#define NN_GRASP_COEFF (int)((double)RAND_MAX/0.8)

enum NNThreadsOption {
    NN_ST,
    NN_MT
};

Solution NearestNeighbor(Instance *inst);

Solution NearestNeighborGrasp(Instance *inst, double timeLimit, enum NNThreadsOption option);

#endif // TSP_HEURISTICS_H