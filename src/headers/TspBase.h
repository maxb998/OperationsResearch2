#ifndef TSP_BASE
#define TSP_BASE

#include <stdlib.h>

#define MAX_THREADS 32
#define USE_APPROXIMATED_DISTANCES 1

#define LOG_LEVEL LOG_LVL_DEBUG

// size of avx vector. 4 is vector of doubles 64bits, 8 is vector of floats 32bits
#define AVX_VEC_SIZE 8

#define SMALLX 1e-6
#define EPSILON 1e-9

enum edgeWeightType{
	EUC_2D, // euclidean distance 2d
	MAN_2D, // manhattan distance 2d
	MAX_2D, // maximum distance 2d
	CEIL_2D, // euclidean 2d rounded up
	ATT, // special distance for problems att48 and att532
	EUC_3D, // euclidean distance 3d
	MAN_3D, // manhattan distance 3d
	MAX_3D, // maximum distance 3d
    GEO, // geographical distance
	XRAY1, // special distance for crystallography problems v1
	XRAY2, // special distance for crystallography problems v2
	EXPLICIT, // weights are specified in the file
	SPECIAL // special type of distance documented elsewhere
};

// data structures
typedef struct
{
    int edgeWeightType;
    int randomSeed;			// if no value has been specified as argument GRASP won't be enabled
	int roundWeights;
    char inputFile[1000];
	char name[200];
    int nThreads;	// if no value has been specified as argument its default value is the number of processors in the machine
} Parameters;

typedef struct
{
    // data
    size_t nNodes;
    float *X;
    float *Y;

	// edge cost mat data
	float *edgeCostMat;

    Parameters params;
} Instance;

typedef struct
{
    float bestCost;    // best solution found cost

	double execTime;

	// the solution
	float *X;
	float *Y;
	// Stores in order of visit along the path, the original indexes (the ones contained in the Instance pointed by instance) of the points
	int *indexPath;

	Instance *instance;
} Solution;

#endif // TSP_BASE