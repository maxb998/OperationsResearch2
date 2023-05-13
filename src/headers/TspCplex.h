#ifndef TSP_CPLEX
#define TSP_CPLEX

#include "TspBase.h"
#include "EdgeCostFunctions.h"
#include <cpxconst.h> // contains only basic data (avoids user to indirectly include cplex.h as a whole)

typedef struct
{
	CPXENVptr env;
	CPXLPptr lp;
	Instance *inst;
} CplexData;

// Contains useful data for the callback. It will go in the parameter cbhandle of the callback call.
typedef struct
{
	Solution *sol;
	CplexData *cpx;

	// Number of columns of the solution
	int ncols;
} CallbackData;


//###################################################################################################################################
// TSP_CPLEX_BASE
//###################################################################################################################################

CplexData initCplexData(Instance *inst);

void destroyCplexData(CplexData * cpxData);

void cplexError(CplexData *cpxData, Instance *inst, Solution *sol, char *line, ...);

size_t xpos(size_t i, size_t j, size_t n);

void cvtCPXtoSuccessors(double *xstar, int ncols, size_t nNodes, int *successors, int *subtoursMap, int *subtourCount);

void cvtSuccessorsToSolution(int *successors, Solution *sol);

// Flag "isBenders" to know what method to use to add the SEC
void setSEC(double *coeffs, int *indexes, CplexData *cpx, int *successors, int *subtoursMap, int subtourCount, int iterNum, Instance *inst, int nCols, int isBenders);

//###################################################################################################################################
// BENDERS
//###################################################################################################################################

Solution benders(Instance *inst, double tlim);

//###################################################################################################################################
// LAZY_CALLBACK
//###################################################################################################################################

// Method that computes the solution using SEC internally
Solution lazyCallback(Instance *inst);

#endif // TSP_CPLEX
