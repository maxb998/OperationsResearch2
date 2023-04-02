#ifndef TSP_IOUTILS
#define TSP_IOUTILS

#include "TspBase.h"

// Read file with .tsp extension according to tsplib specifications, complete with file sintax error checking
void readFile (Instance *inst);

void saveSolution(Solution *sol);

#endif // TSP_IOUTILS