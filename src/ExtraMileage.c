#include "tsp.h"

typedef struct 
{
    Instance *d;
    pthread_mutex_t mutex;
    int *isUncovered; // when an element is set to -1 it means the point with that id is uncovered, if it's 0 than the element is covered
} ThreadInstance;



// Thread task to compute tsp solution using Extra-Mileage procedure using a correctly initialized instance
static void * solveExtraMileageThread(void * arg);


// Finds the hull of the set of points and set it as starting point for Extra-Mileage to work
static inline void quickhull(Instance *d);


static void * quickhullThread();


/* Finds the uppermost, lowermost, leftmost and rightmost points in the set and add them as hull.
* Quickhull algorithm will run starting with this hull.
* Returns number of points added to the hull which varies between 2 and 4 (0 and 1 are bad point-set cases) 
*/
static inline size_t quickhullInitialization(float *x, float *y, size_t n, int *isUncovered, float *xHull, float *yHull);


/* For each point in the set check if it is inside the triagle specified vertex and if that is the case set those points as covered points
* Used by QuickHull algorithm. 
*/
static inline void coverPointsInsideTriangle(float *x, float *y, size_t n, int *isUncovered, float x1, float x2, float x3, float y1, float y2, float y3);


// Returns in a vector areas of the triangles in the input vectors
static inline __m256 detTriangAreaVectorized(__m256 x1, __m256 x2, __m256 x3, __m256 y1, __m256 y2, __m256 y3);


// Returns 0 if there are still uncovered nodes, otherwise returns 1
static  int checkHullCoverageComplete(int *isUncovered, size_t n, size_t memSize);






double solveExtraMileage(Instance *d)
{
    // check if data has been initialized
    // TODO

    // allocate memory for solution
    d->solution.bestSolution = aligned_alloc(32, d->nodesCount + sizeof(int));
}

static void * solveExtraMileageThread(void * arg)
{

}

//######################################################################################################################################################################################################
// QuickHull algorithm functions
//######################################################################################################################################################################################################
static inline void quickhull(Instance *d)
{
    size_t n = d->nodesCount, nMem = d->edgeCost.rowSizeMem;

    int *isUncovered = aligned_alloc(32, nMem * sizeof(int));
    memset(isUncovered, -1, nMem * sizeof(int)); // set all elements as uncovered (-1 to avoids confusion in vectorized programming)

    // initialize hull data
    size_t hullSize = 0; // worst case(not interesting since tsp would be obvious) this goes up to nodesCount
    float *xHull = aligned_alloc(32, nMem * sizeof(float));
    float *yHull = aligned_alloc(32, nMem * sizeof(float));

    // get initial part of the hull in some kind of a smart way
    quickhullInitialization(d->X, d->Y, n, isUncovered, xHull, yHull);

    ThreadInstance th = { .d = d, .isUncovered = isUncovered };
    pthread_mutex_init(&th.mutex, NULL);

    // create threads
    pthread_t threads[MAX_THREADS];
    for (size_t i = 0; i < d->params.threadsCount; i++)
    {
        /* code */
    }
    

    // Preparations: COMPLETED -> begin quickhull iterative process of finding distantmost uncovered point and adding it to the hull
    do
    {
        // find distantmost point from hull points
            
    } while (checkHullCoverageComplete(isUncovered, d->nodesCount, d->edgeCost.rowSizeMem) == 0);

    free(isUncovered);
    free(xHull);
    free(yHull);
}

static inline size_t quickhullInitialization(float *x, float *y, size_t n, int *isUncovered, float *xHull, float *yHull)
{
    size_t nMem = n + AVX_VEC_SIZE - n % AVX_VEC_SIZE;

    // to start the algoritm select first the points that have greatest x, greatest y, smallest x, smallest y
    // create variables that will be used
    __m256 xMaxVec = _mm256_setzero_ps(), xMinVec = _mm256_set1_ps( (float)INFINITY ), yMaxVec = _mm256_setzero_ps(), yMinVec = _mm256_set1_ps( (float)INFINITY );
    __m256i xMaxVecID = _mm256_set_epi32(7,6,5,4,3,2,1,0), xMinVecID = _mm256_set_epi32(7,6,5,4,3,2,1,0), yMaxVecID = _mm256_set_epi32(7,6,5,4,3,2,1,0), yMinVecID = _mm256_set_epi32(7,6,5,4,3,2,1,0);
    register __m256i idsVec = _mm256_set_epi32(7,6,5,4,3,2,1,0), ones = _mm256_set1_epi32(1), mask;

    for (size_t i = 0; i < n - AVX_VEC_SIZE; i += AVX_VEC_SIZE)
    {
        // load new x coords
        register __m256 dataVec = _mm256_load_ps(&x[i]);

        // compare get mask to select correct ids to update in xMaxVecID and only than update xMaxVec with new max(if any) and at the end update the xMaxVecID
        mask = _mm256_castps_si256(_mm256_cmp_ps(dataVec, xMaxVec, _CMP_GT_OS));
        xMaxVec = _mm256_max_ps(xMaxVec, dataVec);
        xMaxVecID = _mm256_blendv_epi8(xMaxVecID, idsVec, mask);

        // same procedure as for xMaxVec
        mask = _mm256_castps_si256(_mm256_cmp_ps(dataVec, xMinVec, _CMP_LT_OS));
        xMinVec = _mm256_min_ps(xMinVec, dataVec);
        xMinVecID = _mm256_blendv_epi8(xMinVecID, idsVec, mask);

        // now same procedure for to find max and min vectors and id vectors in y
        dataVec = _mm256_load_ps(&y[i]);

        // yMaxVec first
        mask = _mm256_castps_si256(_mm256_cmp_ps(dataVec, yMaxVec, _CMP_GT_OS));
        yMaxVec = _mm256_max_ps(yMaxVec, dataVec);
        yMaxVecID = _mm256_blendv_epi8(yMaxVecID, idsVec, mask);

        // same procedure as for xMaxVec
        mask = _mm256_castps_si256(_mm256_cmp_ps(dataVec, yMinVec, _CMP_LT_OS));
        xMinVec = _mm256_min_ps(yMinVec, dataVec);
        xMinVecID = _mm256_blendv_epi8(yMinVecID, idsVec, mask);

        // update vector of ids (increment by one)
        idsVec = _mm256_add_epi32(idsVec, ones);
    }
    
    // now we have 4 vectors one containing in some position the max of x, another the min of x, another the max of y and the last the min of y
    // we also have 4 vectors representing the location of eachpoint (id) in its vector
    // we need to extract these values from the vectors
    // doing in a vectorized way can probably be done but since this is a one time only thing it's ok to not use vectorized instructions
    size_t extremePtsIDs[4];
    float *vecStore = aligned_alloc(32, AVX_VEC_SIZE * sizeof(float));

    // xMax first
    _mm256_store_ps(vecStore, xMaxVec);
    float max = -INFINITY;
    for (size_t i = 0; i < AVX_VEC_SIZE; i++)
        if (max < vecStore[i])
        {
            max = vecStore[i];
            extremePtsIDs[0] = _mm256_extract_epi32(xMaxVecID, (int)i);
        }
    
    // xMin second
    _mm256_store_ps(vecStore, xMinVec);
    float min = INFINITY;
    for (size_t i = 0; i < AVX_VEC_SIZE; i++)
        if (min > vecStore[i])
        {
            min = vecStore[i];
            extremePtsIDs[1] = _mm256_extract_epi32(xMinVecID, (int)i);
        }

    // yMaxVec third
    _mm256_store_ps(vecStore, yMaxVec);
    max = -INFINITY;
    for (size_t i = 0; i < AVX_VEC_SIZE; i++)
        if (max < vecStore[i])
        {
            max = vecStore[i];
            extremePtsIDs[2] = _mm256_extract_epi32(yMaxVecID, (int)i);
        }
    
    // yMinVec last
    _mm256_store_ps(vecStore, yMinVec);
    min = INFINITY;
    for (size_t i = 0; i < AVX_VEC_SIZE; i++)
        if (min > vecStore[i])
        {
            min = vecStore[i];
            extremePtsIDs[3] = _mm256_extract_epi32(yMinVecID, (int)i);
        }
    
    free(vecStore);
    
    // now we have the position of the "extreme" points of the hull
    size_t extremePtsCount = 0;
    for (size_t i = 0; i < 4; i++)
    {
        if (isUncovered[extremePtsIDs[i]] == -1)
        {
            isUncovered[extremePtsIDs[i]] = 0; // add point to the hull
            extremePtsCount++; // increment counter
        }
        else if (i < 3) // remove id of point from extremePtsIDs -> step needed to then compute the points that this starting hull already covers
            extremePtsIDs[extremePtsCount] = extremePtsIDs[i + 1];
    }

    // find points that are covered with this starting hull and set them as such
    if (extremePtsCount >= 3)
        coverPointsInsideTriangle(x, y, n, isUncovered, x[extremePtsIDs[0]], x[extremePtsIDs[1]], x[extremePtsIDs[2]], y[extremePtsIDs[0]], y[extremePtsIDs[1]], y[extremePtsIDs[2]]);
    if (extremePtsCount == 4)
        coverPointsInsideTriangle(x, y, n, isUncovered, x[extremePtsIDs[0]], x[extremePtsIDs[1]], x[extremePtsIDs[3]], y[extremePtsIDs[0]], y[extremePtsIDs[1]], y[extremePtsIDs[3]]);
}

static inline void coverPointsInsideTriangle(float *x, float *y, size_t n, int *isUncovered, float x1, float x2, float x3, float y1, float y2, float y3)
{
    // determine if a point is inside the triangle by computing areas of all possible combination of triangles given 4 points (pt1,pt2,pt3,ptToInvestigate)
    // if the sum of all the area of the triangles that have ptToInvestigate as a vertex is equal (not greater) than the area of the main triangle then the point is inside it

    // use the area of the triangle in determinat form
    float mainTriangArea = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

    // 13 vector register specified
    register __m256 x1Vec = _mm256_set1_ps(x1), x2Vec = _mm256_set1_ps(x2), x3Vec = _mm256_set1_ps(x3);
    register __m256 y1Vec = _mm256_set1_ps(y1), y2Vec = _mm256_set1_ps(y2), y3Vec = _mm256_set1_ps(y3);
    register __m256 xPtsVec, yPtsVec, areaSumVec, mainTriangAreaVec = _mm256_set1_ps(mainTriangArea);
    register __m256i vecIDs = _mm256_set_epi32(7,6,5,4,3,2,1,0), ones = _mm256_set1_epi32(1), result = _mm256_setzero_si256();

    for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
    {
        areaSumVec = _mm256_setzero_ps();

        xPtsVec = _mm256_load_ps(&x[i]);
        yPtsVec = _mm256_load_ps(&y[i]);

        areaSumVec = _mm256_add_ps(areaSumVec, detTriangAreaVectorized(x1Vec, x2Vec, xPtsVec, y1Vec, y2Vec, yPtsVec));
        areaSumVec = _mm256_add_ps(areaSumVec, detTriangAreaVectorized(x1Vec, x3Vec, xPtsVec, y1Vec, y3Vec, yPtsVec));
        areaSumVec = _mm256_add_ps(areaSumVec, detTriangAreaVectorized(x2Vec, x3Vec, xPtsVec, y2Vec, y3Vec, yPtsVec));

        // compare the areas and get -1 to position j if point in position j is inside the mainTriangle -> j is covered
        result = _mm256_castps_si256( _mm256_cmp_ps(mainTriangAreaVec, areaSumVec, _CMP_GT_OS) );

        // store flags
        _mm256_store_si256(&isUncovered[i], result);
    }
}

static inline __m256 detTriangAreaVectorized(__m256 x1, __m256 x2, __m256 x3, __m256 y1, __m256 y2, __m256 y3)
{
    return _mm256_sub_ps(_mm256_mul_ps(_mm256_sub_ps(x1,x3), _mm256_sub_ps(y2,y3)),    _mm256_mul_ps(_mm256_sub_ps(x2,x2), _mm256_sub_ps(y1,y3)));
}

static int checkHullCoverageComplete(int *isUncovered, size_t n, size_t memSize)
{
    // first set placeholder elements to be covered since they do not exist
    if (memSize - n != 0)
        memset(&isUncovered[n], 0, memSize - n);

    // check with and condition if all bits of the array are set to 0 using or logical operator
    register __m256i vec = _mm256_setzero_si256();
    for (size_t i = 0; i < n; i += AVX_VEC_SIZE)
        vec == _mm256_or_si256(vec, _mm256_load_si256(&isUncovered[i]));
    
    // if even one element now inside of vec is different than 0 than the hull does not cover all the points
    for (size_t i = 0; i < AVX_VEC_SIZE; i++)
        //if (_mm256_extract_epi32(vec, i) != 0)
            return 0;
    
    return 1;
}



