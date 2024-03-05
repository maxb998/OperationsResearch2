#include "Tsp.h"
#include "TspCplex.h"

#include <cplex.h>
#include <concorde.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


typedef struct 
{
	CallbackData *cbData;
	void *allocatedMemory;
	bool fromConcorde;
}ConcordePassthroughData;


// Function to post the solution to cplex. Vals and indicies are two allocated portions of memory of ncols elements each.
static inline void PostSolution(CallbackData *cbData, double *vals, int *indicies);
// Function called when a candidate solution has been found(look for subtours and add suitable constraints if necessary)
static inline void candidatePartCallback(CallbackData *cbData);
// Function called during lp relaxation. Finds possible subtours using concorde library and add suitable constrains
static inline void relaxationPartCallback(CallbackData *cbData);
// Function used to add user cuts to cplex environment
int applyCplexUsercut(double cutValue, int cutNcount, int *cutMembers, void *arg);


void BranchAndCut(Solution *sol, double timeLimit)
{
	struct timespec timeStruct;
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

	Instance *inst = sol->instance;

	#ifdef DEBUG
		if (inst->params.cplexWarmStart && (!checkSolution(sol)))
			throwError("Branch&Cut: Input solution is not valid");
	#endif

	CplexData cpx = initCplexData(inst);
	int errCode = 0;

	CallbackData cbData = initCallbackData(&cpx, sol);

	if (inst->params.cplexWarmStart)
		if ((errCode = WarmStart(&cpx, cbData.bestSuccessors)) != 0)
			throwError("Branch&Cut: error on WarmStart with code %d", errCode);

	if (CPXsetintparam(cpx.env, CPX_PARAM_MIPCBREDLP, CPX_OFF))
		throwError("Branch&Cut: error on CPXsetinitparam(CPX_PARAM_MIPCBREDLP)");

	CPXLONG contextMask = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (inst->params.cplexUsercuts) contextMask = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;

	if (CPXcallbacksetfunc(cpx.env, cpx.lp, contextMask, genericCallbackCandidate, &cbData))
		throwError("Branch&Cut: error on CPXsetlazyconstraintcallbackfunc");
	
	if (CPXsetintparam(cpx.env, CPX_PARAM_THREADS, inst->params.nThreads))
		throwError("Branch&Cut: error on CPXsetintparam(CPX_PARAM_THREADS)");

	CPXsetdblparam(cpx.env, CPX_PARAM_TILIM, timeLimit);

	if (CPXmipopt(cpx.env, cpx.lp))
		throwError("Branch&Cut: output of CPXmipopt != 0");

	// asses if the best solution has been found based on whether there is time remainig from the time limit or not(probably not good method)
	clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double currentTime = cvtTimespec2Double(timeStruct);
	if (currentTime - startTime < inst->params.tlim)
	{
		LOG(LOG_LVL_NOTICE, "Branch & Cut found the optimal solution");
	}
	
	cvtSuccessorsToSolution(cbData.bestSuccessors, sol);
	//sol->cost = cbData.bestCost;

	destroyCallbackData(&cbData);

	clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
	sol->execTime += cvtTimespec2Double(timeStruct) - startTime;
}

int CPXPUBLIC genericCallbackCandidate(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle ) 
{
	CallbackData *cbData = (CallbackData *) userhandle;
	cbData->context = context;

	if (contextid & CPX_CALLBACKCONTEXT_RELAXATION)
		relaxationPartCallback(cbData);
	else if (contextid & CPX_CALLBACKCONTEXT_CANDIDATE)
		candidatePartCallback(cbData);

	return 0;
}

CallbackData initCallbackData(CplexData *cpx, Solution *sol)
{
	CallbackData cbData = {
		.inst = cpx->inst,
		.ncols = CPXgetnumcols(cpx->env, cpx->lp),
		.iterNum = 0,
		.bestCost = -1LL,
		.elist = NULL
	};

	if (cbData.ncols <= 0)
		throwError("initCallbackData: CPXgetnumcols returned 0. Cpx.env is probably empty");

	cbData.bestSuccessors = malloc((cpx->inst->nNodes + AVX_VEC_SIZE) * sizeof(int));
	if (cbData.bestSuccessors == NULL)
		throwError("initCallbackData: Failed to allocate memory");

	if (cbData.inst->params.cplexUsercuts)
	{
		cbData.elist = malloc(cbData.ncols * 2 * sizeof(int));
		if (!cbData.elist)
			throwError("initCallbackData: Failed to allocate memory");
		
		for (int i = 0, k = 0; i < cbData.inst->nNodes; i++)
		{
			for (int j = i+1; j < cbData.inst->nNodes; j++)
			{
				cbData.elist[k] = i;
				k++;
				cbData.elist[k] = j;
				k++;
			}
		}
	}

	pthread_mutex_init(&cbData.mutex, NULL);

	cbData.bestCost = sol->cost;
	if (cbData.inst->params.cplexWarmStart)
		cvtSolutionToSuccessors(sol, cbData.bestSuccessors);

	return cbData;
}

void destroyCallbackData(CallbackData *cbData)
{
	free(cbData->bestSuccessors);
	free(cbData->elist);
	cbData->bestSuccessors = cbData->elist = NULL;
	pthread_mutex_destroy(&cbData->mutex);
}

static inline void candidatePartCallback(CallbackData *cbData)
{
	// Stores current vector x from CPLEX
	double *xstar = malloc(cbData->ncols * (sizeof(double) + sizeof(int)));
	if (xstar == NULL)
		throwError("Branch&Cut: could not allocate memory for xstar and indicies");
	// Stores the index of the coefficients that we pass to CPLEX 
	int *indicies = (int*)&xstar[cbData->ncols];

	int errCode;
	if ((errCode = CPXcallbackgetcandidatepoint(cbData->context, xstar, 0, cbData->ncols-1, NULL)) != 0) 
		throwError("Branch&Cut: CPXcallbackgetcandidatepoint failed with code %d", errCode);
	
	Instance *inst = cbData->inst;

	SubtoursData sub = initSubtoursData(inst->nNodes);
	
	cvtCPXtoSuccessors(xstar, cbData->ncols, inst->nNodes, &sub);

	pthread_mutex_lock(&cbData->mutex);
	if (inst->params.mode == MODE_BRANCH_CUT) // avoids too much output in hard fixing and local branching
		LOG(LOG_LVL_DEBUG, "Branch & Cut: iteration %d subtours detected %d", cbData->iterNum, sub.subtoursCount);
	cbData->iterNum++;
	pthread_mutex_unlock(&cbData->mutex);

	__uint128_t cost = -1LL;
	if(sub.subtoursCount > 1)
	{
		double * coeffs = xstar;
		int errCode = setSEC(coeffs, indicies, NULL, cbData, &sub, 0, inst, cbData->ncols);
		if (errCode != 0) // ignore this type of errors
			LOG(LOG_LVL_WARNING, "Branch & Cut: setSEC failed with code %d. Continuing...", errCode);

		cost = PatchingHeuristic(&sub, inst);
	}
	else
		cost = computeSuccessorsSolCost(sub.successors, inst);

	// set new incumbent if found
	if (cost < cbData->bestCost)
	{
		pthread_mutex_lock(&cbData->mutex);
		if (cost < cbData->bestCost)
		{
			swapElems(cbData->bestSuccessors, sub.successors)

			LOG(LOG_LVL_LOG, "Branch & Cut: Found new best solution with cost: %lf", cvtCost2Double(cost));

			cbData->bestCost = cost;

			// post solution to cplex
			#ifdef DEBUG
				if (sub.subtoursCount < 1)
					throwError("BranchAndCut callback: The number of subtours in the solution intended for posting is %d while it must be == 1");
			#endif
			if (inst->params.cplexSolPosting)
				PostSolution(cbData, xstar, indicies);
		}
		pthread_mutex_unlock(&cbData->mutex);
	}

	destroySubtoursData(&sub);
	free(xstar);
}

static inline void PostSolution(CallbackData *cbData, double *vals, int *indicies)
{
	#ifdef DEBUG
		if (!checkSuccessorSolution(cbData->inst, cbData->bestSuccessors) != 0)
			throwError("Successor solution to be posted is incorrect");
	#endif

	for (int i = 0; i < cbData->ncols; i++)
		vals[i] = 0.0;
	for (int i = 0; i < cbData->ncols; i++)
		indicies[i] = i;
	for (int i = 0; i < cbData->inst->nNodes; i++)
	{
		int pos = xpos(i, cbData->bestSuccessors[i], cbData->inst->nNodes);
		vals[pos] = 1.0;
	}

	int strat = CPXCALLBACKSOLUTION_NOCHECK;
	#ifdef DEBUG
		strat = CPXCALLBACKSOLUTION_CHECKFEAS;
	#endif

	int errCode = CPXcallbackpostheursoln(cbData->context, cbData->ncols, indicies, vals, cvtCost2Double(cbData->bestCost), strat);
	if (errCode != 0) // ignore this type of errors
		LOG(LOG_LVL_WARNING, "Branch & Cut solution posting failed with code %d. Continuing...", errCode);
	else
		LOG(LOG_LVL_DEBUG, "Branch & Cut: new solution has been posted correctly");
}


static inline void relaxationPartCallback(CallbackData *cbData)
{
	// Stores current vector x from CPLEX
	size_t allocSize = cbData->ncols * (sizeof(double) + sizeof(int));
	double *xstar = malloc(allocSize);
	if (xstar == NULL)
		throwError("Branch&Cut: could not allocate memory for xstar and indicies");
	// Stores the index of the coefficients that we pass to CPLEX 
	//int *indicies = (int*)&xstar[cbData->ncols];

	int errCode;
	if ((errCode = CPXcallbackgetrelaxationpoint(cbData->context, xstar, 0, cbData->ncols-1, NULL)) != 0) 
		throwError("Branch&Cut: CPXcallbackgetrelaxationpoint failed with code %d", errCode);
	
	int ncomp = 0, *compsCount = NULL, *comps = NULL, ecount = cbData->ncols;

	if ((errCode = CCcut_connect_components(cbData->inst->nNodes, ecount, cbData->elist, xstar, &ncomp, &compsCount, &comps)) != 0)
		throwError("Branch&Cut: CCcut_connect_components failed with code %d", errCode);
	
	ConcordePassthroughData passthroughData = { .cbData=cbData, .allocatedMemory=(void*)xstar, .fromConcorde=false };

	if (ncomp == 1)
	{

		LOG(LOG_LVL_EVERYTHING, "Cutting Single component");
		passthroughData.fromConcorde = true;
		if ((errCode = CCcut_violated_cuts(cbData->inst->nNodes, cbData->ncols, cbData->elist, xstar, 1.9, applyCplexUsercut, &passthroughData)) != 0)
			LOG(LOG_LVL_WARNING, "Branch&Cut - relaxationPartCallback: CCcut_violated_cuts returned error code %d", errCode);
	}
	else if (ncomp > 1)
	{
		for (int i = 0, compPos = 0; i < ncomp; i++)
		{
			applyCplexUsercut(0.0, compsCount[i], &comps[compPos], &passthroughData);
			compPos += compsCount[i];
		}
	}


	free(xstar);
	free(compsCount);
	free(comps);
}

int applyCplexUsercut(double cutValue, int cutNcount, int *cutMembers, void *arg)
{
	ConcordePassthroughData *passthroughData = arg; 
	CallbackData *cbData = passthroughData->cbData;
	int n = cbData->inst->nNodes;

	int cutEcount = cutNcount * (cutNcount - 1) / 2;

	double *coeffs = (double*)passthroughData->allocatedMemory;
	int *index = (int*)&coeffs[cutEcount];
	double rhs = (double)(cutNcount - 1);
	char sense = 'L';
	int purgeable = CPX_USECUT_FILTER;
	int matbeg = 0;
	int local = 0;

	for (int i = 0; i < cutEcount; i++)
		coeffs[i] = 1.0;
	
	int k = 0;
	for (int i = 0; i < cutNcount; i++)
		for (int j = i+1; j < cutNcount; j++, k++)
			index[k] = xpos(cutMembers[i], cutMembers[j], n);
	

	int errCode = CPXcallbackaddusercuts(cbData->context, 1, cutEcount, &rhs, &sense, &matbeg, index, coeffs, &purgeable, &local);
	if (errCode != 0) // ignore this type of errors
		LOG(LOG_LVL_DEBUG, "Branch&Cut - applyCplexUsercut: CPXcallbackaddusercuts returned error code %d. Continuing...", errCode);
	
	if (passthroughData->fromConcorde)
		LOG(LOG_LVL_EVERYTHING, "Added CONCORDE usercut");
	else
		LOG(LOG_LVL_EVERYTHING, "Added manual usercut");

	return 0;
}
