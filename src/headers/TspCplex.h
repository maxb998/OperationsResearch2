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

	int *elist;

	__uint128_t bestCost;
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
int xpos(int i, int j, int n);

/*!
* @brief Convert Cplex solution to a SubtoursData format
* @attention subData.successors and subData.subtoursMap must point to allocated memory spaces of size nNodes
* @param xstar Cplex solution array in cplex format
* @param ncols Number of elements composing xstar
* @param inst Pointer to the instance
* @param subData Output of conversion is stored here
*/
void cvtCPXtoSuccessors(double *xstar, int ncols, Instance *inst, SubtoursData *subData);

/*!
* @brief Compute the cost of a solution in "successors" form
* @param successors Pointer to solution in "successors" form
* @param inst Pointer to the instance of the problem
* @result Cost of the solution
*/
__uint128_t computeSuccessorsSolCost(int *successors, Instance *inst);

/*!
* @brief Convert a successors array in a standard Solution struct
* @attention Successors must not contain any subtour
* @param successors Array of successors
* @param cost Cost of the successor solution
* @param sol Pointer to already initialized solution that points to the target Instance
*/
void cvtSuccessorsToSolution(int *successors, __uint128_t cost, Solution *sol);

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
* @param cpx Cplex Environment and Linear Problem Pointers. Must be NULL if in a callback
* @param cbData Callback Data. Must be NULL if not in callback
* @param subData Pointer to Data realtive to subtours
* @param iterNum Integer used for logging information and cplex constraints naming
* @param inst Pointer to Instance of the problem
* @param nCols Number of Columns/Variables defined by the problem
* @result 0 if everything run smoothly, Cplex Error code is returned otherwise
*/
int setSEC(double *coeffs, int *indexes, CplexData *cpx, CPXCALLBACKCONTEXTptr context, SubtoursData *subData, int iterNum, Instance *inst, int nCols);

/*!
* @brief "Fix" a successors-based solution that presents subtours by merging them in the best possible way.
* @param sub Data relative to the subtours
* @param inst Instance of the problem
* @attention Mutex only locks bestSuccessorsSol and bestCost, the other ones must be thread-local data
* @result Cost of the Repaired Solution
*/
__uint128_t PatchingHeuristic(SubtoursData *sub, Instance *inst);

/*!
* @brief Check correctness of subtoursData
* @param inst Instance of the problem
* @param sub subtoursData to check
* @result true if everything is correct, false otherwise
*/
bool checkSubtoursData(Instance *inst, SubtoursData *sub);

/*!
* @brief Check Solution in the form of an array of successors
* @param inst Instance of the problem
* @param successors Successors array
* @result true if everything is correct, false otherwise
*/
bool checkSuccessorSolution(Instance *inst, int *successors);

/*!
* @brief Add feasible solution in successors form to cplex warm start set
* @param cpx Pointer to the initialized cplex data
* @param sol Successor solution to add to the warm start set. Must be feasible
* @result Returns 0 on success, Cplex Error code on failure
*/
int WarmStart(CplexData *cpx, int *successors);

/*!
* @brief Allocate memory and initialize new SubtoursData struct
* @param n Number of elements in solution
* @result Initialized SubtoursData
*/
SubtoursData initSubtoursData(int n);

/*!
* @brief Free all memory allocated by SubtoursData struct
* @param subData Pointer to SubtoursData struct target
*/
void destroySubtoursData(SubtoursData *subData);

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
// LAZY_CALLBACK - BRANCH AND CUT
//###################################################################################################################################

/*!
* @brief Implement Branch and Cut method in Cplex using a generic callback to add SEC to the problem.
* @param sol Pointer to a valid solution that will be used as warm start as well as output
* @param tlim Time limit for benders
*/
void BranchAndCut(Solution *sol, double tlim);

// cplex callback to that solves tsp instance
int CPXPUBLIC genericCallbackCandidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle );

/*!
* @brief Initializes CallbackData by setting values and allocating memory
* @param cpx Pointer to CplexData
* @param sol Pointer to warm start solution
*/
CallbackData initCallbackData(CplexData *cpx, Solution *sol);

/*!
* @brief Free memory allocated inside a CallbackData type and destroy the mutex
* @param cbData Pointer to the CallbackData target
*/
void destroyCallbackData(CallbackData *cbData);

//###################################################################################################################################
// HARD_FIXING
//###################################################################################################################################

/*!
* @brief Run Hard Fixing matheuristic on sol
* @param sol Feasible solution to optimize with Hard-Fixing
* @param timeLimit Time limit
*/
void HardFixing(Solution *sol, double timeLimit);


//###################################################################################################################################
// LOCAL_BRANCHING
//###################################################################################################################################

/*!
* @brief Run Local Branching matheuristic on sol
* @param sol Pointer to feasible solution to optimize
* @param timeLimit Time limit
*/
void LocalBranching(Solution *sol, double timeLimit);

#endif // TSP_CPLEX