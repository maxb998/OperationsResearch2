#ifndef TSP_BASE
#define TSP_BASE

#include <stdlib.h>
#include <stdbool.h>

#define MAX_THREADS 64

// size of avx vector. 4 is vector of doubles 64bits, 8 is vector of floats 32bits
#define AVX_VEC_SIZE 8

#define SMALLX 1e-6
#define EPSILON 1e-9


#define COMPUTE_OPTION_AVX 0
#define COMPUTE_OPTION_BASE 1 
#define COMPUTE_OPTION_USE_COST_MATRIX 2

#define COMPUTATION_TYPE COMPUTE_OPTION_AVX

// Amount of "best" elements to save when using grasp almostbest option and NOT using COMPUTE_OPTION_AVX
#define BASE_GRASP_BEST_SAVE_BUFFER_SIZE 4


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
    EM_INIT_FARTHEST_POINTS
};

struct VnsKickSize
{
	int Min;
	int Max;
};

struct GeneticParams
{
	int populationSize;
	int crossoverAmount;
	int mutationAmount;
	int reintroductionAmount;
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
	double tlim;

	enum GraspOption graspType;
	double graspChance;
	enum NNFirstNodeOptions nnFirstNodeOption;
	enum EMInitType emInitOption;

	enum Mode metaheurInitMode;
	int metaRestartThreshold;
	int tabuTenureSize;
	struct VnsKickSize vnsKickSize;
	struct GeneticParams geneticParams;
	int annealingTemperature;

	enum Mode matheurInitMode;	
	enum HardFixPolicy hardFixPolicy;
    
	bool use2Opt;

    int randomSeed;
	int nThreads; // if no value has been specified as argument its default value is the number of processors in the machine
	bool roundWeights;
	bool showPlot;
	bool saveSolution;
	enum LogLevel logLevel;

	enum EdgeWeightType edgeWeightType;
	char name[200];
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
    __uint128_t cost;    // best solution found cost

	double execTime;

	// Stores in order of visit along the path, the original indexes (the ones contained in the Instance pointed by instance) of the points
	int *indexPath;

	Instance *instance;
} Solution;


#endif // TSP_BASE