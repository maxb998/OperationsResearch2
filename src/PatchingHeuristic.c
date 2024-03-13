#include "TspCplex.h"
#include "Tsp.h"

typedef struct
{
	int index0;
	int index1;
	bool invertOrientation;
} MergingData;


static inline MergingData findBestSubtourMerge(SubtoursData *sub, int subtoursCount, Instance *inst);

__uint128_t PatchingHeuristic(SubtoursData *sub, Instance *inst)
{
	int n = inst->nNodes;

	int subtoursCount = sub->subtoursCount;
	while (subtoursCount > 1)
	{
		MergingData md = findBestSubtourMerge(sub, subtoursCount, inst);

		if (!md.invertOrientation)
			swapElems(sub->successors[md.index0], sub->successors[md.index1])
		else
		{
			int last = sub->successors[md.index1];
			int previous = sub->successors[md.index1];
			int current = sub->successors[previous];
			int next = sub->successors[current];
			do
			{
				sub->successors[current] = previous;
				previous = current;
				current = next;
				next = sub->successors[next];
			} while (current != last);
			
			sub->successors[last] = sub->successors[md.index0];
			sub->successors[md.index0] = md.index1;
		}

		// now adjust subtoursMap
		int sub2Map = sub->subtoursMap[md.index1];
		for (int i = 0; i < n; i++)
		{
			if (sub->subtoursMap[i] == sub2Map)
				sub->subtoursMap[i] = sub->subtoursMap[md.index0];
			else if (sub->subtoursMap[i] > sub2Map)
				sub->subtoursMap[i]--;
		}

		subtoursCount--;
	}

	__uint128_t cost = computeSuccessorsSolCost(sub->successors, inst);

	#ifdef DEBUG
		if (!checkSuccessorSolution(inst, sub->successors))
			throwError("Patching Heurisitc: Produced successors is not valid");
	#endif

	return cost;
}

static inline MergingData findBestSubtourMerge(SubtoursData *sub, int subtoursCount, Instance *inst)
{
	#if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
		float *X = inst->X, *Y =inst->Y;
	#endif

	float min = INFINITY;

	MergingData retVal = { .index0=0, .index1=0, .invertOrientation=false };

	for (int subtourID = 0; subtourID < subtoursCount-1; subtourID++)
	{
		int first1 = 0;
			while (sub->subtoursMap[first1] != subtourID)
				first1++;

		for (int subtourToCompare = subtourID+1; subtourToCompare < subtoursCount; subtourToCompare++)
		{
			int first2 = 0;
			while (sub->subtoursMap[first2] != subtourToCompare)
				first2++;
			
			bool finish1 = false;
			int i = first1;
			while ((!finish1) || (i != sub->successors[first1]))
			{
				if (sub->successors[i] == first1) finish1 = true;
				
				// successor of i'th node
				int succI = sub->successors[i];

				bool finish2 = false;
				int j = first2;
				while ((!finish2) || (j != sub->successors[first2]))
				{
					if (sub->successors[j] == first2) finish2 = true;

					// successor of j'th node
					int succJ = sub->successors[j];
					
					#if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
						float cost = computeEdgeCost(X[i], Y[i], X[succJ], Y[succJ], inst) + computeEdgeCost(X[succI], Y[succI], X[j], Y[j], inst);
					#elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
						float cost = inst->edgeCostMat[i * (size_t)inst->nNodes + succJ] + inst->edgeCostMat[succI * (size_t)inst->nNodes + j];
					#endif

					if (cost < min)
					{
						min = cost;
						retVal.index0 = i; retVal.index1 = j;
						retVal.invertOrientation = false;
					}
					
					#if ((COMPUTATION_TYPE == COMPUTE_OPTION_AVX) || (COMPUTATION_TYPE == COMPUTE_OPTION_BASE))
						cost = computeEdgeCost(X[i], Y[i], X[j], Y[j], inst) + computeEdgeCost(X[succI], Y[succI], X[succJ], Y[succJ], inst);
					#elif (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
						cost = inst->edgeCostMat[i * (size_t)inst->nNodes + j] + inst->edgeCostMat[succI * (size_t)inst->nNodes + succJ];
					#endif

					if (cost < min)
					{
						min = cost;
						retVal.index0 = i; retVal.index1 = j;
						retVal.invertOrientation = true;
					}

					j = succJ;
				}
				i = succI;
			}
		}
	}

	return retVal;
}

