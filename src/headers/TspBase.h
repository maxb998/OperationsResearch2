#ifndef TSP_BASE
#define TSP_BASE

#include <stdlib.h>
#include <stdbool.h>

#define MAX_THREADS 64
#define USE_APPROXIMATED_DISTANCES 1

// size of avx vector. 4 is vector of doubles 64bits, 8 is vector of floats 32bits
#define AVX_VEC_SIZE 8

#define SMALLX 1e-6
#define EPSILON 1e-9

enum LogLevel{
	LOG_LVL_ERROR,
	LOG_LVL_CRITICAL,
	LOG_LVL_WARNING,
	LOG_LVL_NOTICE,
	LOG_LVL_LOG,
	LOG_LVL_DEBUG,
	LOG_LVL_EVERYTHING
};

enum EdgeWeightType{
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

enum Mode{
	MODE_NONE=-1,
	MODE_NN,
	MODE_EM,
	MODE_TABU,
	MODE_VNS,
	MODE_ANNEALING,
	MODE_GENETIC,
	MODE_BENDERS,
	MODE_BRANCH_CUT,
	MODE_HARDFIX,
	MODE_LOCAL_BRANCHING
};

enum GraspOption{
	GRASP_NONE=-1,
	GRASP_ALMOSTBEST,
	GRASP_RANDOM
};

enum NNFirstNodeOptions{
	NN_FIRST_RANDOM,
    NN_FIRST_TRYALL
};

enum EMInitType {
    EM_INIT_RANDOM,
    EM_INIT_FARTHEST_POINTS,
    EM_INIT_HULL // won't be working for now
};

enum EMOptions {
    EM_OPTION_AVX,
    EM_OPTION_BASE,
    EM_OPTION_USE_COST_MATRIX
};

enum _2OptOptions
{
    _2OPT_AVX_ST,
    _2OPT_BASE,
    _2OPT_PRECOMPUTED_COSTS,
    _2OPT_AVX_MT
};

enum HardFixPolicy {
	HARDFIX_POLICY_RANDOM,
	HARDFIX_POLICY_SMALLEST
};

// data structures
typedef struct
{
	char inputFile[1000];
	enum Mode mode;

	enum GraspOption graspType;
	bool use2Opt;
	double tlim;

	enum NNFirstNodeOptions nnFirstNodeOption;
	enum EMInitType emInitOption;

	enum Mode metaheurInitMode;
	enum Mode matheurInitMode;
	enum Mode warmStartMode;

	enum HardFixPolicy hardFixPolicy;
    
    int randomSeed;
	int nThreads; // if no value has been specified as argument its default value is the number of processors in the machine
	bool roundWeights;
	bool showPlot;
	bool saveSolution;
	enum LogLevel logLevel;

	enum EdgeWeightType edgeWeightType;
	char name[200];
	double graspChance;
} Parameters;

typedef struct
{
    // data
    int nNodes;
    float *X;
    float *Y;

	// edge cost mat data
	float *edgeCostMat;

    Parameters params;
} Instance;

typedef struct
{
    double cost;    // best solution found cost

	double execTime;

	// the solution
	float *X;
	float *Y;
	// Stores in order of visit along the path, the original indexes (the ones contained in the Instance pointed by instance) of the points
	int *indexPath;

	Instance *instance;
} Solution;


#endif // TSP_BASE