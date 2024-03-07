#ifndef TSP_BASE
#define TSP_BASE

#include <stdlib.h>
#include <stdbool.h>

#define MAX_THREADS 64

// size of avx vector. 4 is vector of doubles 64bits, 8 is vector of floats 32bits
#define AVX_VEC_SIZE 8

#define SMALLX 1e-6
#define EPSILON 1e-9

#define USE_REDUCED_DISTANCE_SET // only use EUC_2D, CEIL_2D and ATT distances since those are the only ones in the selected 81 tsplib .tsp files

#define COMPUTE_OPTION_AVX 0
#define COMPUTE_OPTION_BASE 1 
#define COMPUTE_OPTION_USE_COST_MATRIX 2

#define COMPUTATION_TYPE COMPUTE_OPTION_AVX

// #define DEBUG

// Amount of "best" elements to save when using grasp almostbest option and NOT using COMPUTE_OPTION_AVX
#define BASE_GRASP_BEST_SAVE_BUFFER_SIZE 4


enum LogLevel{
	LOG_LVL_FATAL,
	LOG_LVL_ERROR,
	LOG_LVL_WARN,
	LOG_LVL_NOTICE,
	LOG_LVL_INFO,
	LOG_LVL_DEBUG,
	LOG_LVL_TRACE
};


enum EdgeWeightType{
	#ifndef USE_REDUCED_DISTANCE_SET
	MAN_2D, // manhattan distance 2d
	MAX_2D, // maximum distance 2d
	#endif
	EUC_2D, // euclidean distance 2d
	CEIL_2D, // euclidean 2d rounded up
	ATT // special distance for problems att48 and att532
};

enum Mode{
	MODE_NONE=0,
	MODE_NN=1,
	MODE_EM=2,
	MODE_TABU=4,
	MODE_VNS=8,
	MODE_ANNEALING=16,
	MODE_GENETIC=32,
	MODE_BENDERS=64,
	MODE_BRANCH_CUT=128,
	MODE_HARDFIX=256,
	MODE_LOCAL_BRANCHING=512
};

enum GraspOption{
	GRASP_NONE=0,
	GRASP_ALMOSTBEST=1,
	GRASP_RANDOM=2
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
	bool cplexPatching;
	bool cplexWarmStart;
	bool cplexSolPosting;
	bool cplexUsercuts;
	enum HardFixPolicy hardFixPolicy;
    
	bool use2Opt;
	bool use3Opt;

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

	#if (COMPUTATION_TYPE == COMPUTE_OPTION_USE_COST_MATRIX)
	float *edgeCostMat;
	#endif

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