#ifndef EDGE_COST_FUNCTIONS
#define EDGE_COST_FUNCTIONS

#include "TspBase.h"
#include <math.h>
#include <immintrin.h>


static inline float noSquaredRootEdgeCost (float x1, float y1, float x2, float y2, Instance *inst)
{
	register float costSquared;
	if(inst->params.edgeWeightType == MAN_2D)
		costSquared = fabsf(x1-x2) + fabsf(y1-y2);
	else if (inst->params.edgeWeightType == MAX_2D)
	{
		costSquared = fabsf(x1 - x2);
		register float diffY = fabsf(y1 - y2);
		if (diffY > costSquared)
			costSquared = diffY;
	}
	else // either euc2d or att
	{
		register float diffX = x1 - x2, diffY = y1 - y2;
		costSquared = (diffX * diffX) + (diffY * diffY);
		if (inst->params.edgeWeightType == ATT)
			costSquared /= 10.F;
	}

	return costSquared;
}


static inline float computeEdgeCost (float x1, float y1, float x2, float y2, Instance *inst)
{
	register float cost = noSquaredRootEdgeCost(x1, y1, x2, y2, inst);

	if ((inst->params.edgeWeightType == EUC_2D) || (inst->params.edgeWeightType == ATT))
		cost = sqrtf(cost);

	if (inst->params.roundWeights)
	{
		if (inst->params.edgeWeightType == ATT)
			cost = ceilf(cost);
		else
			cost = roundf(cost);
	}
	
	return cost;
}



// ###########################################################################################################

static inline __m256 noSquaredRootEdgeCost_VEC (__m256 x1, __m256 y1, __m256 x2, __m256 y2, Instance *inst)
{
	register __m256 costSquared;
	if(inst->params.edgeWeightType == MAN_2D)
	{
		register __m256 xDiff, yDiff, fabsfMask = _mm256_set1_ps(-0.0F);

		xDiff = _mm256_sub_ps(x1, x2);
		yDiff = _mm256_sub_ps(y1, y2);
		xDiff = _mm256_andnot_ps(fabsfMask, xDiff);
		yDiff = _mm256_andnot_ps(fabsfMask, yDiff);

		costSquared = _mm256_add_ps(xDiff, yDiff);
	}
	else if (inst->params.edgeWeightType == MAX_2D)
	{
		register __m256 xDiff, yDiff, fabsfMask = _mm256_set1_ps(-0.0F);

		xDiff = _mm256_sub_ps(x1, x2);
		yDiff = _mm256_sub_ps(y1, y2);
		xDiff = _mm256_andnot_ps(fabsfMask, xDiff);
		yDiff = _mm256_andnot_ps(fabsfMask, yDiff);

		costSquared = _mm256_max_ps(xDiff, yDiff);
	}
	else // either euc2d or att
	{
		register __m256 xDiff = _mm256_sub_ps(x1, x2);
		register __m256 yDiff = _mm256_sub_ps(y1, y2);
		costSquared = _mm256_add_ps( _mm256_mul_ps(xDiff, xDiff), _mm256_mul_ps(yDiff, yDiff));

		if (inst->params.edgeWeightType == ATT)
			costSquared = _mm256_mul_ps(costSquared, _mm256_set1_ps(0.1F));
	}

	return costSquared;
}

static inline __m256 computeEdgeCost_VEC (__m256 x1, __m256 y1,  __m256 x2, __m256 y2, Instance *inst)
{
	register __m256 costVec = noSquaredRootEdgeCost_VEC (x1, y1, x2, y2, inst);

	if ((inst->params.edgeWeightType == EUC_2D) || (inst->params.edgeWeightType == ATT))
		costVec = _mm256_sqrt_ps(costVec);

	if (inst->params.roundWeights)
	{
		if (inst->params.edgeWeightType == ATT)
			return _mm256_ceil_ps(costVec);
		else
			return _mm256_round_ps(costVec, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC); // round to nearest, and suppress exceptions
	}
	
	return costVec;
}

static inline __m256 computeEdgeCost_VEC_APPROX(__m256 x1, __m256 y1,  __m256 x2, __m256 y2, Instance *inst)
{
	register __m256 costVec = noSquaredRootEdgeCost_VEC (x1, y1, x2, y2, inst);

	if ((inst->params.edgeWeightType == EUC_2D) || (inst->params.edgeWeightType == ATT))
		costVec = _mm256_rcp_ps(_mm256_rsqrt_ps(costVec));

	if (inst->params.roundWeights)
	{
		if (inst->params.edgeWeightType == ATT)
			return _mm256_ceil_ps(costVec);
		else
			return _mm256_round_ps(costVec, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC); // round to nearest, and suppress exceptions
	}
	
	return costVec;
}

#endif // EDGE_COST_FUNCTIONS
