#ifndef TSP_CPLEX
#define TSP_CPLEX

#include "TspBase.h"
#include "EdgeCostFunctions.h"
#include <cpxconst.h> // contains only basic data (avoids user to indirectly include cplex.h as a whole)
#include <pthread.h>

typedef struct
{
	CPXENVptr env;
	CPXLPptr lp;
	Instance *inst;
} CplexData;

typedef struct 
{
	int *successors;
	int *subtoursMap;
	int subtoursCount;
} SubtoursData;

typedef struct
{
	Instance *inst;

	int ncols; // Number of columns of the matrix inside cplex linear problem
	int iterNum;
	pthread_mutex_t mutex;

	double bestCost;
	int *bestSuccessors;
} CallbackData;

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
* @brief Compute the cost of a solution in "successors" form
* @param successors Pointer to solution in "successors" form
* @param inst Pointer to the instance of the problem
* @result Cost of the solution
*/
double computeSuccessorsSolCost(int *successors, Instance *inst);

/*!
* @brief Convert a successors array in a standard Solution struct
* @attention Successors must not contain any subtour
* @param successors Array of successors
* @param sol Pointer to already initialized solution that points to the target Instance
*/
void cvtSuccessorsToSolution(int *successors, Solution *sol);

/*!
* @brief Convert a Solution to a successors array
* @param sol Solution to convert
* @param successors Output memory
*/
void cvtSolutionToSuccessors(Solution *sol, int* successors);

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
* @result 0 if everything run smoothly, Cplex Error code is returned otherwise
*/
int setSEC(double *coeffs, int *indexes, CplexData *cpx, CPXCALLBACKCONTEXTptr context, SubtoursData *subData, int iterNum, Instance *inst, int nCols, bool isBenders);

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
bool checkSuccessorSolution(Instance *inst, int *successors);

/*!
* @brief Add feasible solution in successors form to cplex warm start set
* @param cpx Pointer to the initialized cplex data
* @param sol Successor solution to add to the warm start set. Must be feasible
* @result Returns 0 on success, Cplex Error code on failure
*/
int WarmStart(CplexData *cpx, int *successors);

//###################################################################################################################################
// BENDERS
//###################################################################################################################################

/*!
* @brief Apply benders method to the tsp instance inst. Uses Repair Heurstic between calls of CPXmipopt to provide a solution if time limit does not allow to find optimal solution.
* @param sol Pointer to a valid solution that will be used as warm start as well as output
* @param tlim Time limit for benders
*/
void benders(Solution *sol, double tlim);

//###################################################################################################################################
// LAZY_CALLBACK
//###################################################################################################################################

/*!
* @brief Implement Branch and Cut method in Cplex using a generic callback to add SEC to the problem.
* @param sol Pointer to a valid solution that will be used as warm start as well as output
* @param tlim Time limit for benders
*/
void BranchAndCut(Solution *sol, double tlim);

// cplex callback to that solves tsp instance
int CPXPUBLIC genericCallbackCandidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle );

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
