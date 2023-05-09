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


//###################################################################################################################################
// TSP_CPLEX_BASE
//###################################################################################################################################

CplexData initCplexData(Instance *inst);

void destroyCplexData(CplexData * cpxData);

size_t xpos(size_t i, size_t j, size_t n);

//###################################################################################################################################
// BLENDERS
//###################################################################################################################################

Solution blenders(Instance *inst, double tlim);


#endif // TSP_CPLEX