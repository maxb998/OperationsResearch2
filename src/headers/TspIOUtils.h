#ifndef TSP_IOUTILS
#define TSP_IOUTILS

#include "TspBase.h"

// Read file with .tsp extension according to tsplib specifications, complete with file sintax error checking.
// Returns time elapsed while reading file
double readFile (Instance *inst);

void saveSolution(Solution *sol);

/*Plot solution using gnuplot. Does NOT check for errors on input
 * d	-> Instance to plot
 * plotPixelSize	-> string: Plot window size in pixel specified with format: "<WIDTH>,<HEIGHT>"
 * pointColor -> string: Color of the circle representing the point, eg "black" or "red"
 * tourPointColor -> string: Color of the 'X' on top of the point circle of color pointColor. Format and types is the same as for pointColor
 * pointSize -> int: Size of the points
 * printIndex -> int: set to 1 to print index of each point as label on plot
*/
void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize, const int printIndex);

#endif // TSP_IOUTILS