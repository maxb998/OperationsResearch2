#include "tsp.h"

typedef struct
{
    Instance *d;
    pthread_mutex_t nodeLock;
    pthread_mutex_t saveLock;
    int startingNode;
}ThreadedInstance;

double NearestNeighbour(Instance *d, int configuration);

static void * computeSingleNN(void *thInst);

int findSuccessor(Instance *d, int *uncoveredNodes, int node, float *pathCost);

double NearestNeighbour(Instance *d, int configuration)
{
    // first we check the number of processors to now how many threads we want to create
    FILE *commandPipe;
    char *command = "nproc";
    char temp[10];
    commandPipe = (FILE*)popen(command, "r");
    fgets(temp, sizeof(temp), commandPipe);
    int numProcessors = atoi(temp);

    // we create and initialize the threaded instance
    ThreadedInstance thInst = {.d = d, .startingNode = 0};
    pthread_mutex_init(&thInst.nodeLock, NULL);
    pthread_mutex_init(&thInst.saveLock, NULL);
    
    pthread_t threads[numProcessors];
    for(int i = 0; i < numProcessors; i++)
    {
        pthread_create(&threads[i], NULL, computeSingleNN, &thInst);
        LOG(LOG_LVL_LOG, "Nearest Neighbour : Thread %d CREATED", i);
    }

    //double threadsTime[numProcessors];    // SEGMENTATION FAULT BELOW
    double *returnPtr;
    for(int i = 0; i < numProcessors; i++)
    {
        pthread_join(threads[i], (void **)&returnPtr);
        //threadsTime[i] = *returnPtr;      // THROWS SEGMENTATION FAULT
        free(returnPtr);

        LOG(LOG_LVL_LOG, "Thread %d finished in x time", i);
    }




    pthread_mutex_destroy(&thInst.nodeLock);
    pthread_mutex_destroy(&thInst.saveLock);
    return -1;
}

static void * computeSingleNN(void *thInst)
{
    ThreadedInstance *th = (ThreadedInstance *)thInst;

    // we create an array that keeps track of the successors along the path
    int * iterationPath = malloc(th->d->nodesCount * sizeof(int));
    // we create an array that indicates if a node has alredy been visited
    int * uncoveredNodes = malloc(th->d->nodesCount * sizeof(int));

    // we want the threads to repeat the computation for every node
    // for this we use a mutex on startingNode until it reaches nodesCount
    while((pthread_mutex_lock(&th->nodeLock) == 0) && (th->startingNode < th->d->nodesCount))
    {
        int node = th->startingNode;
        th->startingNode++;
        // after incrementing the value of startingNode we can unlock the mutex
        pthread_mutex_unlock(&th->nodeLock);

        // set all elements of uncoveredNodes to zero
        memset(uncoveredNodes, 0, th->d->nodesCount * sizeof(int));

        // initialize the cost of the path to zero
        float pathCost = 0;

        // we set the starting node as visited
        uncoveredNodes[node] = 1;
        int currentNode = node;
        int successor;
        for(int i = 0; i < th->d->nodesCount-1; i++)    // for n nodes we want to run this loop n-1 times, at the end we set as successor of the last node the starting node
        {
            // find the successor node
            // pathCost is also updated in this method
            successor = findSuccessor(th->d, uncoveredNodes, currentNode, &pathCost);
            // set the successor in the path
            iterationPath[currentNode] = successor;
            // update current node
            currentNode = successor;
        }
        // at the end we set the starting node as successor of the last one to close the circuit
        iterationPath[currentNode] = node;

        // to check if we have to update the best solution we use another mutex
        if((pthread_mutex_lock(&th->saveLock) == 0) && (th->d->solution.bestCost > pathCost))
        {
            th->d->solution.bestCost = pathCost;
            int * temp = th->d->solution.bestSolution;
            th->d->solution.bestSolution = iterationPath;
            iterationPath = temp;
        }
        pthread_mutex_unlock(&th->saveLock);
    }
    pthread_mutex_unlock(&th->nodeLock);
    return 0;
}

int findSuccessor(Instance *d, int *uncoveredNodes, int node, float *pathCost)
{
    float bestDistance = INFINITY;
    int currentBestNode = -1;
    // check if we are working with rounded weights
    if(d->params.roundWeights == 0)
    {
        for(int i = 0; i < d->nodesCount; i++)
        {
        // check if the node has alredy been visited, if it's different from the current node
        // and if the distance is better than the best seen
            if(uncoveredNodes[i] == 0 && i != node && d->edgeCost.mat[(d->edgeCost.rowSizeMem)*node + i] < bestDistance)
            {
                currentBestNode = i;
                bestDistance = d->edgeCost.mat[(d->edgeCost.rowSizeMem)*node + i];
                uncoveredNodes[i] = 1;
            }
        }
    }else
    {
        for(int i = 0; i < d->nodesCount; i++)
        {
        // here we are working with rounded weights, which are stored in d.edgeCost.roundedMat
            if(uncoveredNodes[i] == 0 && i != node && d->edgeCost.roundedMat[(d->edgeCost.rowSizeMem)*node + i] < bestDistance)
            {
                currentBestNode = i;
                bestDistance = d->edgeCost.mat[(d->edgeCost.rowSizeMem)*node + i];
                uncoveredNodes[i] = 1;
            }
        }
    }
    *pathCost += bestDistance;
    return currentBestNode;
}