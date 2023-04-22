#include "TspCplex.h"
#include "EdgeCostFunctions.h"

#include <cplex.h>

typedef struct
{
	CPXENVptr env;
	CPXLPptr lp;
	Instance *inst;
} CplexData;


int xpos(size_t i, size_t j, size_t n);

static CplexData initCplexData(Instance *inst);

static void destroyCplexData(CplexData * cpxData);


Solution blenders(Instance *inst, double tLimSec)
{
    Solution sol = newSolution(inst);

    CplexData cpx = initCplexData(inst);

	

	destroyCplexData(&cpx);

    return sol;
}

int xpos(size_t i, size_t j, size_t n)
{ 
	if ( i == j )
		throwError(NULL, NULL,"xpos: i == j");
	if ( i > j )
	{
		register size_t temp;
		swapElems(i,j,temp);
	}

	int pos = i * n + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}

static CplexData initCplexData(Instance *inst)
{
	CplexData cpxData;

	int errno = 0;
	cpxData.env = CPXopenCPLEX(&errno);
	if (errno)
		throwError(inst, NULL, "buildCPXModel: error at CPXopenCPLEX with code %d", errno);
	cpxData.lp = CPXcreateprob(cpxData.env, &errno, "TSP");
	if (errno)
		throwError(inst, NULL, "buildCPXModel: error at CPXcreateprob with code %d", errno);

    size_t n = inst->nNodes;
    enum edgeWeightType ewt = inst->params.edgeWeightType;
    int roundFlag = inst->params.roundWeightsFlag;

	char binary = CPX_BINARY; 

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

    // add binary var.s x(i,j) for i < j  
	
	for ( int i = 0; i < n; i++ )
	{
		for ( int j = i+1; j < n; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = computeEdgeCost(inst->X[i], inst->Y[i], inst->X[j], inst->Y[j], ewt, roundFlag); // cost == distance
			double ub = 1.0;
			if ( CPXnewcols(cpxData.env, cpxData.lp, 1, &obj, NULL, &ub, &binary, cname) ) 
				throwError(inst, NULL, "initCplexData: wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(cpxData.env, cpxData.lp)-1 != xpos(i,j,n) )
				throwError(inst, NULL, "initCplexData: wrong position for x var.s");
		}
	} 

    // add the degree constraints 

	int *index = (int *) calloc(n, sizeof(int));
	double *value = (double *) calloc(n, sizeof(double));

	for ( int h = 0; h < n; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < n; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h,n);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(cpxData.env, cpxData.lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) 
			throwError(inst, NULL, "initCplexData: CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	return cpxData;
}

static void destroyCplexData(CplexData * cpxData)
{
	CPXfreeprob(cpxData->env, &cpxData->lp);
    CPXcloseCPLEX(&cpxData->env);
}
