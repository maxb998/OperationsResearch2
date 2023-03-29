#include "tsp.h"

typedef struct 
{
    Instance inst;
    pthread_mutex_t mutex;
    int *isUncovered; // when an element is set to -1 it means the point with that id is uncovered, if it's 0 than the element is covered
} ThreadInstance;
