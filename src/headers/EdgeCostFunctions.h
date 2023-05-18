#ifndef EDGE_COST_FUNCTIONS
#define EDGE_COST_FUNCTIONS

#include "TspBase.h"
#include <math.h>
#include <immintrin.h>

enum costComputationType
{
	FAST_SQUARED,		// use squared weight when possible (avoids square root function which takes a lot of computing power)
	MODERATE_APPROX,	// use square root approximation (approximation error of less than 0.1% (max error at around 0.073%))
	SLOW_EXACT			// use standard square root function (basically accurate)
};

#define COMPUTATION_TYPE SLOW_EXACT



static inline float euclideanCostSquared2D(float x1, float y1, float x2, float y2)
{
    register float diffX = x1 - x2, diffY = y1 - y2;
    return (diffX * diffX) + (diffY * diffY);
}

static inline float euclideanCost2D(float x1, float y1, float x2, float y2)
{
    return sqrtf(euclideanCostSquared2D(x1, y1, x2, y2));
}

static inline float manhattanCost2D(float x1, float y1, float x2, float y2)
{
    return (fabsf(x1-x2) + fabsf(y1-y2));
}

static inline float maximumCost2D(float x1, float y1, float x2, float y2)
{
    register float diffX = x1 - x2, diffY = y1 - y2;
    if (diffX > diffY)
        return diffX;
    else
        return diffY;
}

static inline float attCostSquared2D(float x1, float y1, float x2, float y2)
{
    return (euclideanCostSquared2D(x1, y1, x2, y2) / 10.0F);
}

static inline float attCost2D(float x1, float y1, float x2, float y2)
{
    return sqrtf(attCostSquared2D(x1, y1, x2, y2));
}

static inline float squaredEdgeCost (float x1, float y1, float x2, float y2, enum EdgeWeightType edgeWgtType)
{
	switch (edgeWgtType)
	{
	case EUC_2D:
		return euclideanCostSquared2D(x1, y1, x2, y2);
		break;
	
	case MAN_2D:
		return manhattanCost2D(x1, y1, x2, y2);
		break;

	case MAX_2D:
		return maximumCost2D(x1, y1, x2, y2);
		break;

	case ATT:
		return attCostSquared2D(x1, y1, x2, y2);
		break;
	
	default:  // euclidean 2D if cost is not known
		return euclideanCostSquared2D(x1, y1, x2, y2);
		break;
	}
}

static inline float exactEdgeCost (float x1, float y1, float x2, float y2, enum EdgeWeightType edgeWgtType)
{
	switch (edgeWgtType)
	{
	case EUC_2D:
		return euclideanCost2D(x1, y1, x2, y2);
		break;
	
	case MAN_2D:
		return manhattanCost2D(x1, y1, x2, y2);
		break;

	case MAX_2D:
		return maximumCost2D(x1, y1, x2, y2);
		break;

	case ATT:
		return attCost2D(x1, y1, x2, y2);
		break;
	
	default:  // euclidean 2D if cost is not known
		return euclideanCost2D(x1, y1, x2, y2);
		break;
	}
}

static inline float roundEdgeCost(float edgeCost, enum EdgeWeightType edgeWgtType)
{
	if (edgeWgtType == ATT)
		return ceilf(edgeCost);
	else
		return floorf(edgeCost);
}

static inline float computeEdgeCost (float x1, float y1, float x2, float y2, enum EdgeWeightType edgeWgtType, int roundFlag)
{
	register float cost;
	switch (COMPUTATION_TYPE)
	{
	case FAST_SQUARED:
		cost = squaredEdgeCost (x1, y1, x2, y2, edgeWgtType);
		break;
	
	case MODERATE_APPROX:
		cost = exactEdgeCost (x1, y1, x2, y2, edgeWgtType);
		break;

	case SLOW_EXACT:
		cost = exactEdgeCost (x1, y1, x2, y2, edgeWgtType);
		break;
	}

	if (roundFlag == 1)
		cost = roundEdgeCost (cost, edgeWgtType);
	
	return cost;
}





// Return vector containing the euclidean cost squared. Fastest
static inline __m256 euclideanCostSquared2D_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 xDiff, yDiff, dist;

    xDiff = _mm256_sub_ps(x1, x2);
    yDiff = _mm256_sub_ps(y1, y2);
    dist = _mm256_add_ps(_mm256_mul_ps(xDiff, xDiff), _mm256_mul_ps(yDiff, yDiff));

    return dist;
}
// Return vector containing the most accurate euclidean cost. Slow
static inline __m256 euclideanCost2D_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    return _mm256_sqrt_ps(euclideanCostSquared2D_VEC(x1,y1,x2,y2));
}
// Return vector containig the approximations of euclidean costs. Maximum relative error is 3*2^-12. Fast euclidean cost
static inline __m256 euclideanCost2DFastApprox_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    return _mm256_rcp_ps(_mm256_rsqrt_ps(euclideanCostSquared2D_VEC(x1,y1,x2,y2)));
}

// https://stackoverflow.com/questions/63599391/find-fabsfolute-in-avx
static inline __m256 manhattanCost2D_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 xDiff, yDiff, dist, fabsfMask = _mm256_set1_ps(-0.0F);
    
    xDiff = _mm256_sub_ps(x1, x2);
    yDiff = _mm256_sub_ps(y1, y2);
    xDiff = _mm256_andnot_ps(fabsfMask, xDiff);
    yDiff = _mm256_andnot_ps(fabsfMask, yDiff);
    dist = _mm256_add_ps(xDiff, yDiff);

    return dist;
}

static inline __m256 maximumCost2D_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 xDiff, yDiff, dist;

    xDiff = _mm256_sub_ps(x1, x2);
    yDiff = _mm256_sub_ps(y1, y2);
    dist = _mm256_max_ps(xDiff, yDiff);

    return dist;
}

static inline __m256 attCostSquared2D_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    register __m256 dist, vec10 = _mm256_set1_ps(10.F);
    dist = euclideanCostSquared2D_VEC(x1, y1, x2, y2);

    return _mm256_div_ps(dist, vec10);
}

static inline __m256 attCost2D_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    return _mm256_sqrt_ps(attCostSquared2D_VEC(x1,y1,x2,y2));
}

static inline __m256 attCost2DFastApprox_VEC(__m256 x1, __m256 y1, __m256 x2, __m256 y2)
{
    return _mm256_rcp_ps(_mm256_rsqrt_ps(attCostSquared2D_VEC(x1,y1,x2,y2)));
}

// ###########################################################################################################

static inline __m256 squaredEdgeCost_VEC (__m256 x1, __m256 y1, __m256 x2, __m256 y2, enum EdgeWeightType edgeWgtType)
{
	switch (edgeWgtType)
	{
	case EUC_2D:
		return euclideanCostSquared2D_VEC(x1, y1, x2, y2);
		break;
	
	case MAN_2D:
		return manhattanCost2D_VEC(x1, y1, x2, y2);
		break;

	case MAX_2D:
		return maximumCost2D_VEC(x1, y1, x2, y2);
		break;

	case ATT:
		return attCostSquared2D_VEC(x1, y1, x2, y2);
		break;
	
	default:  // euclidean 2D if cost is not known
		return euclideanCostSquared2D_VEC(x1, y1, x2, y2);
		break;
	}
}

static inline __m256 exactEdgeCost_VEC (__m256 x1, __m256 y1,  __m256 x2, __m256 y2, enum EdgeWeightType edgeWgtType)
{
	switch (edgeWgtType)
	{
	case EUC_2D:
		return euclideanCost2D_VEC(x1, y1, x2, y2);
		break;
	
	case MAN_2D:
		return manhattanCost2D_VEC(x1, y1, x2, y2);
		break;

	case MAX_2D:
		return maximumCost2D_VEC(x1, y1, x2, y2);
		break;

	case ATT:
		return attCost2D_VEC(x1, y1, x2, y2);
		break;
	
	default:  // euclidean 2D if cost is not known
		return euclideanCost2D_VEC(x1, y1, x2, y2);
		break;
	}
}

static inline __m256 fastEdgeCost_VEC (__m256 x1, __m256 y1,  __m256 x2, __m256 y2, enum EdgeWeightType edgeWgtType)
{
	switch (edgeWgtType)
	{
	case EUC_2D:
		return euclideanCost2DFastApprox_VEC(x1, y1, x2, y2);
		break;
	
	case MAN_2D:
		return manhattanCost2D_VEC(x1, y1, x2, y2);
		break;

	case MAX_2D:
		return maximumCost2D_VEC(x1, y1, x2, y2);
		break;

	case ATT:
		return attCost2DFastApprox_VEC(x1, y1, x2, y2);
		break;
	
	default:  // euclidean 2D if cost is not known
		return euclideanCost2DFastApprox_VEC(x1, y1, x2, y2);
		break;
	}
}

static inline __m256 roundEdgeCost_VEC (__m256 costs, enum EdgeWeightType edgeWgtType)
{
	if (edgeWgtType == ATT)
		return _mm256_ceil_ps(costs);
	else
		return _mm256_floor_ps(costs);
}

static inline __m256 computeEdgeCost_VEC (__m256 x1, __m256 y1,  __m256 x2, __m256 y2, enum EdgeWeightType edgeWgtType, int roundFlag)
{
	register __m256 costVec;
	switch (COMPUTATION_TYPE)
	{
	case FAST_SQUARED:
		costVec = squaredEdgeCost_VEC (x1, y1, x2, y2, edgeWgtType);
		break;
	
	case MODERATE_APPROX:
		costVec = fastEdgeCost_VEC (x1, y1, x2, y2, edgeWgtType);
		break;
	
	case SLOW_EXACT:
		costVec = exactEdgeCost_VEC (x1, y1, x2, y2, edgeWgtType);
		break;
	}

	if (roundFlag == 1)
		costVec = roundEdgeCost_VEC(costVec, edgeWgtType);
	
	return costVec;
}

#endif // EDGE_COST_FUNCTIONS
