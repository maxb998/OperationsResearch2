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

typedef struct 
{
	int *successors;
	int *subtoursMap;
	int subtoursCount;
} SubtoursData;


//###################################################################################################################################
// TSP_CPLEX_BASE
//###################################################################################################################################

/*!
* @brief Generate CplexData struct using a tsp instance as input
* @param inst Instance to initialize within cplex
* @result Initialized cplex basic structure
*/
CplexData initCplexData(Instance *inst);

/*!
* @brief Destroys cplex enviroment and instanced problem. (Free memory)
* @param cpxData Pointer to CplexData.
*/
void destroyCplexData(CplexData * cpxData);

/*!
* @brief Print error message and frees any memory from each pointer that is not NULL. Then exits the program with code 1.
* @param cpxData Pointer to CplexData that will be destroyed before exiting
* @param inst Pointer to Instance that will be destroyed before exiting.
* @param sol Pointer to Solution that will be destroyed before exiting.
* @param line Error message and format
*/
void cplexError(CplexData *cpxData, Instance *inst, Solution *sol, char *line, ...);

/*!
* @brief Print error message and returns 1. Useful when getting errors inside a callback.
* @param line Error message and format
* @result 1
*/
int callbackError(char *line, ...);

/*!
* @brief Convert a simple edge coordinates (given by 2 nodes indexes) to an index for the cplex edges array format. 
* @attention i and j must be different
* @param i First node of edge
* @param j Second node of edge
* @param n Number of nodes
* @result Index of edge selected by i and j in cplex format
*/
size_t xpos(size_t i, size_t j, size_t n);

/*!
* @brief Convert Cplex solution to a SubtoursData format
* @attention subData.successors and subData.subtoursMap must point to allocated memory spaces of size nNodes
* @param xstar Cplex solution array in cplex format
* @param ncols Number of elements composing xstar
* @param nNodes Number of nodes in the instance
* @param subData Output of conversion is stored here
*/
void cvtCPXtoSuccessors(double *xstar, int ncols, size_t nNodes, SubtoursData *subData);

/*!
* @brief Convert a successors array in a standard Solution struct
* @attention Successors must not contain any subtour
* @param successors Array of successors
* @param sol Pointer to already initialized solution that points to the target Instance
*/
void cvtSuccessorsToSolution(int *successors, Solution *sol);

/*!
* @brief Set subtour elimination constraints given a solution that have more than one subtour.
* @param coeffs Allocated memory of nCols doubles used to add SEC
* @param indexes Allocated memory of nCols integers used to add SEC
* @param cpx Cplex Environment and Linear Problem Pointers. Can be NULL if in a callback
* @param context Contex of Callback if inside a Callback. Can be NULL if not in callback
* @param subData Pointer to Data realtive to subtours
* @param iterNum Integer used for logging information and cplex constraints naming
* @param inst Pointer to Instance of the problem
* @param nCols Number of Columns/Variables defined by the problem
* @param isBenders Flag to decide which function to use to add constraints to the problem
*/
void setSEC(double *coeffs, int *indexes, CplexData *cpx, CPXCALLBACKCONTEXTptr context, SubtoursData *subData, int iterNum, Instance *inst, int nCols, int isBenders);

/*!
* @brief "Fix" a successors-based solution that presents subtours by merging them in the best possible way.
* @param sub Data relative to the subtours
* @param inst Instance of the problem
* @attention Mutex only locks bestSuccessorsSol and bestCost, the other ones must be thread-local data
* @result Cost of the Repaired Solution
*/
double PatchingHeuristic(SubtoursData *sub, Instance *inst);

/*!
* @brief Check Solution in the form of an array of successors
* @param inst Instance of the problem
* @param successors Successors array
* @result 0 if everything is correct. 1 if successors is not correct
*/
int checkSuccessorSolution(Instance *inst, int *successors);

//###################################################################################################################################
// BENDERS
//###################################################################################################################################

/*!
* @brief Apply benders method to the tsp instance inst. Uses Repair Heurstic between calls of CPXmipopt to provide a solution if time limit does not allow to find optimal solution.
* @param inst Instance of the tsp problem to solve
* @param tlim Time limit for benders
* @result Optimal solution or solution built using repair heuristic if time limit wasn't long enough
*/
Solution benders(Instance *inst, double tlim);

//###################################################################################################################################
// LAZY_CALLBACK
//###################################################################################################################################

/*!
* @brief Implement Branch and Cut method in Cplex using a generic callback to add SEC to the problem.
* @param inst Instance of tsp to solve
* @param tlim Time limit for benders
* @result Solution of the problem obtained by branch and cut method
*/
Solution BranchAndCut(Instance *inst, double tlim);


//###################################################################################################################################
// HARD_FIXING
//###################################################################################################################################

/*!
* @brief 
* @param sol Feasible solution to optimize with Hard-Fixing
* @param tlim Time limit
*/
void HardFixing(Solution *sol, double fixingAmount, enum HardFixPolicy policy, double tlim);

#endif // TSP_CPLEX
