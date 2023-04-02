#include "ExtraMileage.h"
#include "EdgeCostFunctions.h"
#include "TspUtilities.h"

#include "pthread.h"
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time


typedef struct 
{
    Instance inst;
    pthread_mutex_t mutex;
    int *isUncovered; // when an element is set to -1 it means the point with that id is uncovered, if it's 0 than the element is covered
} ThreadInstance;
