#include "tsp.h"

typedef struct
{
    Instance *inst;
    pthread_mutex_t nodeLock;
    pthread_mutex_t saveLock;
    int startingNode;
}ThreadedInstance;

static void * threadNN(void *thInst);

/* Returns:
    1. Index of closest unvisited node 
    2. -1 in case of error
    3. -2 as flag if we alredy visited all the nodes
*/
static inline int findSuccessor(Instance *inst, int *uncoveredNodes, int node);

// Given two nodes node1 and node2 updates the cost of the path.
static inline void updatePathCost(Instance *inst, double *pathCost, int node1, int node2);

double NearestNeighbour(Instance *inst)
{
    // we initialize the seed if it has been passed as argument
    if(inst->params.randomSeed != -1) srand(inst->params.randomSeed);
    // we create and initialize the threaded instance
    ThreadedInstance thInst = {.inst = inst, .startingNode = 0};
    pthread_mutex_init(&thInst.nodeLock, NULL);
    pthread_mutex_init(&thInst.saveLock, NULL);
    
    pthread_t threads[inst->params.threadsCount];
    for(int i = 0; i < inst->params.threadsCount; i++)
    {
        pthread_create(&threads[i], NULL, threadNN, &thInst);
        LOG(LOG_LVL_LOG, "Nearest Neighbour : Thread %d CREATED", i);
    }

    //double threadsTime[numProcessors];    // SEGMENTATION FAULT BELOW
    double *returnPtr;
    for(int i = 0; i < inst->params.threadsCount; i++)
    {
        pthread_join(threads[i], (void **)&returnPtr);
        //threadsTime[i] = *returnPtr;      // THROWS SEGMENTATION FAULT
        free(returnPtr);

        LOG(LOG_LVL_LOG, "Nearest Neighbour : Thread %d finished in x time", i);
    }




    pthread_mutex_destroy(&thInst.nodeLock);
    pthread_mutex_destroy(&thInst.saveLock);
    return 0;
}

static void * threadNN(void *thInst)
{
    ThreadedInstance *th = (ThreadedInstance *)thInst;

    // We create an array that stores the nodes of the solution in order
    // It must contain nodesCount + 1 elements, since the first and the last node are the same
    int * iterationPath = malloc((th->inst->nodesCount+1) * sizeof(int));
    // we create an array that indicates if a node has alredy been visited
    int * uncoveredNodes = malloc(th->inst->nodesCount * sizeof(int));

    // We want the threads to repeat the computation for every node
    // For this we use a mutex on startingNode until it reaches nodesCount
    while((pthread_mutex_lock(&th->nodeLock) == 0) && (th->startingNode < th->inst->nodesCount))
    {
        int node = th->startingNode;
        th->startingNode++;
        // after incrementing the value of startingNode we can unlock the mutex
        pthread_mutex_unlock(&th->nodeLock);

        // set all elements of uncoveredNodes to zero
        memset(uncoveredNodes, 0, th->inst->nodesCount * sizeof(int));
        // DEBUG: reset iterationPath to all zeros
        memset(iterationPath, 0, (th->inst->nodesCount+1) * sizeof(int));

        // initialize the cost of the path to zero
        double pathCost = 0;

        // we set the starting node as visited
        uncoveredNodes[node] = 1;
        // we set the starting node in the path
        int currentNode = iterationPath[0] = iterationPath[th->inst->nodesCount] = node;
        int successor;
        for(int i = 1; i < th->inst->nodesCount; i++)    // for n nodes we want to run this loop n-1 times
        {
            successor = findSuccessor(th->inst, uncoveredNodes, currentNode);


            // Control on validity of successor
            if(successor < 0 || successor >= (int)th->inst->nodesCount)
                throwError(th->inst, "threadNN: error computing successor %d, value returned: %d", i, successor);
            
            // set the successor in the path
            iterationPath[i] = successor;
            // update current node
            currentNode = successor;

            // Update the cost of the path
            updatePathCost(th->inst, &pathCost, iterationPath[i-1], iterationPath[i]);
        }
        updatePathCost(th->inst, &pathCost, iterationPath[th->inst->nodesCount-1], iterationPath[th->inst->nodesCount]);

        // to check if we have to update the best solution we use another mutex
        if((pthread_mutex_lock(&th->saveLock) == 0) && (th->inst->solution.bestCost > pathCost))
        {
            th->inst->solution.bestCost = pathCost;
            int * temp = th->inst->solution.bestSolution;
            th->inst->solution.bestSolution = iterationPath;
            iterationPath = temp;
            if(LOG_LEVEL > LOG_LVL_DEBUG) solutionCheck(th->inst);
            LOG(LOG_LVL_LOG, "Found better solution starting from node %d, cost: %lf", node, pathCost);
        }
        pthread_mutex_unlock(&th->saveLock);
    }
    pthread_mutex_unlock(&th->nodeLock);
    return 0;
}

static inline int findSuccessor(Instance *inst, int *uncoveredNodes, int node)
{
    // to keep track of the closest node
    float bestDistance = INFINITY;
    int bestNode = -1;
    // to keep track of the second closest node (in case that we use GRASP)
    float secondBestDistance = INFINITY;
    int secondBestNode = -1;

    int currentDistance;    // stores the distance of the node that we are checking
 
    for(int i = 0; i < inst->nodesCount; i++)
    {
    // check if the node has alredy been visited, if it's different from the current node
        if(uncoveredNodes[i] == 0 && i != node)
        {
            currentDistance = inst->edgeCost.mat[(inst->edgeCost.rowSizeMem)*node + i];
            // check if the distance is better than the best seen
            if(currentDistance < bestDistance)
            {
                bestNode = i;
                bestDistance = currentDistance;
            // if it's not we check if it's better than the second best distance
            }else if(currentDistance < secondBestDistance)
            {
                secondBestNode = i;
                secondBestDistance = currentDistance;
            }
        }
    }
    // We choose what node of the two best we return if GRASP has been required
    if(inst->params.randomSeed != -1 && rand() > GRASP_COEFF && secondBestNode != -1)
    {
        uncoveredNodes[secondBestNode] = 1;
        return secondBestNode;
    }else
    {    uncoveredNodes[bestNode] = 1;
        return bestNode;
    }
}

static inline void updatePathCost(Instance *inst, double *pathCost, int node1, int node2)
{
    *pathCost += inst->edgeCost.mat[inst->edgeCost.rowSizeMem*node1 + node2];
}