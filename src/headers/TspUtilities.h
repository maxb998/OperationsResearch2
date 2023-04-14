#ifndef UTILITIES
#define UTILITIES

#include "TspBase.h"

void checkSolution(Solution *sol);

float computeSolutionCostVectorizedFloat(Solution *sol);

double computeSolutionCostVectorizedDouble(Solution *sol);

double computeSolutionCost(Solution *sol);

#endif // UTILITIES