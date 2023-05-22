#include "TspCplex.h"
#include "Tsp.h"

static inline void findBestSubtourMerge(SubtoursData *sub, int subtoursCount, Instance *inst, int mergeIndexes[2], int *invertOrientation);

static double computeSuccessorsSolCost(int *successors, Instance *inst);


double PatchingHeuristic(SubtoursData *sub, Instance *inst)
{
	size_t n = inst->nNodes;

	int subtoursCount = sub->subtoursCount;
	while (subtoursCount > 1)
	{
		int invert = 0;
		int mergeIndexes[2] = {0};

		findBestSubtourMerge(sub, subtoursCount, inst, mergeIndexes, &invert);

		if (invert == 0)
		{
			register int temp = sub->successors[mergeIndexes[0]];
			sub->successors[mergeIndexes[0]] = sub->successors[mergeIndexes[1]];
			sub->successors[mergeIndexes[1]] = temp;
		}
		else
		{
			int last = sub->successors[mergeIndexes[1]];
			int previous = sub->successors[mergeIndexes[1]];
			int current = sub->successors[previous];
			int next = sub->successors[current];
			do
			{
				sub->successors[current] = previous;
				previous = current;
				current = next;
				next = sub->successors[next];
			} while (current != last);
			
			sub->successors[last] = sub->successors[mergeIndexes[0]];
			sub->successors[mergeIndexes[0]] = mergeIndexes[1];
		}

		// now adjust subtoursMap
		int sub2Map = sub->subtoursMap[mergeIndexes[1]];
		for (size_t i = 0; i < n; i++)
		{
			if (sub->subtoursMap[i] == sub2Map)
				sub->subtoursMap[i] = sub->subtoursMap[mergeIndexes[0]];
			else if (sub->subtoursMap[i] > sub2Map)
				sub->subtoursMap[i]--;
		}

		subtoursCount--;
	}

	double cost = computeSuccessorsSolCost(sub->successors, inst);

	return cost;
}

static inline void findBestSubtourMerge(SubtoursData *sub, int subtoursCount, Instance *inst, int mergeIndexes[2], int *invertOrientation)
{
	float *X = inst->X, *Y =inst->Y;
	enum EdgeWeightType ewt = inst->params.edgeWeightType ;
	int roundFlag = inst->params.roundWeightsFlag;

	float min = INFINITY;

	for (size_t subtourID = 0; subtourID < subtoursCount-1; subtourID++)
	{
		size_t first1 = 0;
			while (sub->subtoursMap[first1] != subtourID)
				first1++;

		for (size_t subtourToCompare = subtourID+1; subtourToCompare < subtoursCount; subtourToCompare++)
		{
			size_t first2 = 0;
			while (sub->subtoursMap[first2] != subtourToCompare)
				first2++;
			
			int finish1 = 0;
			int i = first1;
			while ((!finish1) || (i != sub->successors[first1]))
			{
				if (sub->successors[i] == first1) finish1 = -1;
				
				// successor of i'th node
				int succI = sub->successors[i];

				int finish2 = 0;
				int j = first2;
				while ((!finish2) || (j != sub->successors[first2]))
				{
					if (sub->successors[j] == first2) finish2 = -1;

					// successor of j'th node
					int succJ = sub->successors[j];

					float cost = computeEdgeCost(X[i], Y[i], X[succJ], Y[succJ], ewt, roundFlag) + computeEdgeCost(X[succI], Y[succI], X[j], Y[j], ewt, roundFlag);

					if (cost < min)
					{
						min = cost;
						mergeIndexes[0] = i; mergeIndexes[1] = j;
						*invertOrientation = 0;
					}
					
					cost = computeEdgeCost(X[i], Y[i], X[j], Y[j], ewt, roundFlag) + computeEdgeCost(X[succI], Y[succI], X[succJ], Y[succJ], ewt, roundFlag);

					if (cost < min)
					{
						min = cost;
						mergeIndexes[0] = i; mergeIndexes[1] = j;
						*invertOrientation = 1;
					}

					j = succJ;
				}
				i = succI;
			}
		}
	}
}

static double computeSuccessorsSolCost(int *successors, Instance *inst)
{
	int n = (int)inst->nNodes;
	enum EdgeWeightType ewt = inst->params.edgeWeightType ;
	int roundFlag = inst->params.roundWeightsFlag;

	double cost = 0.0;

	int i = 0;
	int counter = 0;
	do
	{
		int succ = successors[i];
		cost += computeEdgeCost(inst->X[i], inst->Y[i], inst->X[succ], inst->Y[succ], ewt, roundFlag);
		i = succ;

		if (counter > n)
			throwError(inst, NULL, "computeSuccessorsSolCost: There are subtours inside the successor array even after repair heuristic");
		counter++;
	} while (i != 0);
	
	return cost;
}

