#include "tsp.h"

static enum logLevel globLVL = LOG_LVL_EVERYTHING;

const char * logLevelString [] = {
	"\033[1;31mERR \033[0m", // 0
	"\033[1;35mCRIT\033[0m", // 1
	"\033[1;33mWARN\033[0m", // 2
	"\033[1;36mNOTI\033[0m", // 3
	"\033[1;34mLOG \033[0m", // 4
	"\033[1;32mDEBG\033[0m", // 5
	"\033[1;90mALL \033[0m"	// 6
};

// keywords that we need to find inside .tsp file
const char * keywords[] = {
    "NAME",
    "TYPE",
    "DIMENSION",
    "EDGE_WEIGHT_TYPE",
    "NODE_COORD_SECTION"
};
#define KEYWORDS_COUNT 5

#define NO_KEYWORD_FOUND_ID -1
enum keywordID {
    NAME_KEYWORD_ID,
    TYPE_KEYWORD_ID,
    DIMENSION_KEYWORD_ID,
    EDGE_WEIGHT_TYPE_KEYWORD_ID,
    NODE_COORD_SECTION_KEYWORD_ID
};

// 
const char * wgtTypeStr[] = {
	"EUC_2D", // euclidean distance 2d
	"MAN_2D", // manhattan distance 2d
	"MAX_2D", // maximum distance 2d
	"CEIL_2D", // euclidean 2d rounded up
	"ATT", // special distance for problems att48 and att532
	"EUC_3D", // euclidean distance 3d
	"MAN_3D", // manhattan distance 3d
	"MAX_3D", // maximum distance 3d
    "GEO", // geographical distance
	"XRAY1", // special distance for crystallography problems v1
	"XRAY2", // special distance for crystallography problems v2
	"EXPLICIT", // weights are specified in the file
	"SPECIAL" // special type of distance documented elsewhere
};
#define EDGE_WEIGHT_TYPES_COUNT 12

// file parsing functions

// gets the string specified in the "NAME" keyword in the file
static void getNameFromFile(char * line, int lineSize, char out[], Instance *d);
// checks that the file type is "TSP"
static void checkFileType(char * line, int lineSize, Instance *d);
// get the value related with the keyword "DIMENSION" and returns it as a size_t
static size_t getDimensionFromLine(char * line, int lineSize, Instance *d);
// check that the string associated with "EDGE_WEIGHT_TYPE" is correct and return it as a number
static size_t getEdgeWeightTypeFromLine(char * line, int lineSize, Instance *d);
// Returns the number of processors of the machine
static inline int nProcessors();


void initInstance(Instance *d)
{
    // set pointers to NULL
    d->X = d->Y = d->edgeCost.mat = NULL;
    d->edgeCost.roundedMat = d->solution.bestSolution = NULL;

    // initialize variables to standard value
    d->params.edgeWeightType = -1; // so that it gives error
    d->params.roundWeights = 0;
    d->params.threadsCount = nProcessors();
    d->nodesCount = 0;
    d->solution.bestCost = INFINITY;
    d->params.randomSeed = -1;
    // make strings empty to check for errors (should only need to change first char to 0 but it's nicer this way)
    memset(d->params.inputFile, 0, 1000);
    memset(d->params.name, 0, 200);
}

void freeInstance(Instance *d)
{
    // points
    free(d->X);
    free(d->Y);

    // cost matrix
    free(d->edgeCost.mat);
    free(d->edgeCost.roundedMat);

    // solution
    free(d->solution.bestSolution);
}

int LOG (enum logLevel lvl, char * line, ...)
{
    // check log level
    if (lvl > globLVL) return 0;

    // print log level
    printf("[%s] ", logLevelString[lvl]);

    // print passed message and values
    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);
    
    // add new line at the end
    if (line[strlen(line)-1] != '\n')
        printf("\n");

    return 0;
}

void throwError (Instance *d, char * line, ...)
{
    printf("[%s] ", logLevelString[0]);

    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);

    printf("\n");

    freeInstance(d);

    exit(EXIT_FAILURE);
}

void parseArgs (Instance *d, int argc, char *argv[])
{
    static int roundWeightsFlag = 0;
    static struct option options[] = 
        {
            {"seed", required_argument, 0, 's'},
            {"s", required_argument, 0, 's'},
            {"file", required_argument, 0, 'f'},
            {"f", required_argument, 0, 'f'},
            {"threads", required_argument, 0, 't'},
            {"t", required_argument, 0, 't'},
            {"roundWeigths", no_argument, &roundWeightsFlag, 1},
            {0, 0, 0, 0}
        };
    // set flags in d
    d->params.roundWeights = roundWeightsFlag;


    int option_index = 0;
    int opt;

    while ((opt = getopt_long(argc, argv, "s:f:t:", options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 's':
            d->params.randomSeed = (int)strtol(optarg, NULL, 10);
            break;

        case 'f':
            if (access(optarg, R_OK) != 0)
                throwError(d, "ERROR: File \"%s\" not found\n", optarg);

            strncpy(d->params.inputFile, optarg, strlen(optarg)+1);
            break;
        
        case 't':
            d->params.threadsCount = strtoul(optarg, NULL, 10);
            break;
        
        default:
            abort();
        }
    }

    // check necessary arguments were passed
    if (d->params.inputFile[0] == 0)
        throwError(d, "A file path must be specified with \"-f\" or \"--file\" options");

    LOG(LOG_LVL_NOTICE, "Received arguments:");
    LOG(LOG_LVL_NOTICE,"    Random Seed  = %d", d->params.randomSeed);
    LOG(LOG_LVL_NOTICE,"    Filename     = %s", d->params.inputFile);
    LOG(LOG_LVL_NOTICE,"    Thread Count = %d", d->params.threadsCount);
}


void readFile (Instance *d)
{
    FILE *f = fopen(d->params.inputFile, "r");

    // check if was able to open file
    if (!f) { fprintf (stderr, "failed to open file for reading\n"); exit(EXIT_FAILURE); }

    char *line = NULL;
    size_t lineMemSize;
    int lineSize;
    int keywordsFound[KEYWORDS_COUNT] = { 0 };    // starting value assumes bad file

    size_t keywordsLinesCount = 0; // count the number of lines occupied by the keywords
    int keywordFinished = 0;
    while ((keywordFinished == 0) && ((lineSize = getline(&line, &lineMemSize, f)) != EOF))
    {
        keywordsLinesCount++;
        LOG(LOG_LVL_EVERYTHING, line);

        // check in the first part of each line if we find a useful keyword
        int keywordID = NO_KEYWORD_FOUND_ID;
        for (size_t i = 0; i < KEYWORDS_COUNT; i++)
        {
            size_t keywordLength = strlen(keywords[i]);

            if (strncmp(line, keywords[i], keywordLength) == 0)
            {
                // if inside here means that keyword i has been found
                keywordID = i;
                break;
            }
        }

        if (keywordsFound[keywordID] == 1)
            throwError(d, "Keyword \"%s\" is present more than one time. Check the .tsp file", keywords[keywordID]);

        switch (keywordID)
        {
        case NAME_KEYWORD_ID:
            getNameFromFile(line, lineSize, d->params.name, d);
            // set flag
            keywordsFound[NAME_KEYWORD_ID] = 1;
            break;

        case TYPE_KEYWORD_ID:
            checkFileType(line, lineSize, d);
            // if this point is reached the has a correct type -> set the flag
            keywordsFound[TYPE_KEYWORD_ID] = 1;
            break;
        
        case DIMENSION_KEYWORD_ID:
            d->nodesCount = getDimensionFromLine(line, lineSize, d);
            // set the flag
            keywordsFound[DIMENSION_KEYWORD_ID] = 1;
            break;

        case EDGE_WEIGHT_TYPE_KEYWORD_ID:
            d->params.edgeWeightType = getEdgeWeightTypeFromLine(line, lineSize, d);
            // if this point is reached a correct number for dimension has been taken from the file -> set the flag
            keywordsFound[EDGE_WEIGHT_TYPE_KEYWORD_ID] = 1;
            break;
            
        case NODE_COORD_SECTION_KEYWORD_ID:
            // start reading coordinates of nodes from the next line -> exit the loop
            keywordFinished = 1;
            keywordsFound[NODE_COORD_SECTION_KEYWORD_ID] = 1;
            break;

        case NO_KEYWORD_FOUND_ID:
            break;

        default:
            throwError(d, "Wierd error upon reading keywords from file. Check the code");
            break;
        }
    }

    // check all important keywords have been found
    for (size_t i = 0; i < KEYWORDS_COUNT; i++)
        if (keywordsFound[i] != 1)  // keyword has not been found in the loop above
            throwError(d, "Important keyword \"%s\" has not been found/detected in the tsp file. Check the .tsp file", keywords[i]);

    // allocate memory
    size_t memElemsToAlloc = d->nodesCount + AVX_VEC_SIZE; // allocate more than needed to avoid errors when reading with avx
    d->X = aligned_alloc(32, memElemsToAlloc * sizeof(float));
    d->Y = aligned_alloc(32, memElemsToAlloc * sizeof(float));

    d->solution.bestSolution = aligned_alloc(32, memElemsToAlloc * sizeof(int));

    // fill the memory with data
    size_t i = 0;
    while (((lineSize = getline(&line, &lineMemSize, f)) != EOF) && (strncmp(line, "EOF", 3) != 0))
    {
        // check if the number of data lines are not more than the number specified in dimension
        if (i >= d->nodesCount)
        {
            LOG(LOG_LVL_WARNING, "There may be more elements than the number specified with keyword \"DIMENSION\" in the .tsp file. Ignoring extra data");
            break;
        }

        // get the line number, convert it and check it's ok
        char * endPtr = NULL;
        size_t lineNumber = strtoul(line, &endPtr, 10);
        if (lineNumber == 0)
            throwError(d, "Conversion at line %lu of the line number has gone wrong(must not be 0)", keywordsLinesCount + i + 1);
        
        char * separatorPtr = strchr(line, ' ');
        if (!separatorPtr)
            throwError(d, "Space separator at line %lu of file has not been found");
        if (endPtr != separatorPtr)
            throwError(d, "Conversion at line %lu of the line number has gone wrong(must be an integer separated from other value by a space at the end)", keywordsLinesCount + i + 1);
        
        if (lineNumber != i+1)
            throwError(d, "Line number at line %lu of the file is not what is supposed to be (%lu instead of %lu)", keywordsLinesCount + i + 1, lineNumber, i+1);

        // now we can convert the actual coordinates
        char * xCoordStrPtr = separatorPtr + 1;
        separatorPtr = strchr(xCoordStrPtr, ' ');
        if (!separatorPtr)
            throwError(d, "Space separator at line %lu of file has not been found");

        d->X[i] = strtof(xCoordStrPtr, &endPtr);

        if (xCoordStrPtr == endPtr)
            throwError(d, "Conversion of X coordinate at line %lu has gone wrong. Check the .tsp file", keywordsLinesCount + i + 1);
        if (endPtr != separatorPtr)
            throwError(d, "Conversion of X coordinate at line %lu has gone wrong. There are unwanted characters at the end ", keywordsLinesCount + i + 1);
        if (fabsf(d->X[i]) == HUGE_VALF)
            throwError(d, "Coordinate X at line %lu has caused overflow", keywordsLinesCount + i + 1);

        char * yCoordStrPtr = separatorPtr + 1;
        separatorPtr = strchr(yCoordStrPtr, '\n');
        if (!separatorPtr)
            throwError(d, "'\n' at line %lu of file has not been found");
        
        d->Y[i] = strtof(yCoordStrPtr, &endPtr);

        if (yCoordStrPtr == endPtr)
            throwError(d, "Conversion of Y coordinate at line %lu has gone wrong. Check the .tsp file", keywordsLinesCount + i + 1);
        if (endPtr != separatorPtr)
            throwError(d, "Conversion of Y coordinate at line %lu has gone wrong. There are unwanted characters at the end ", keywordsLinesCount + i + 1);
        if (fabsf(d->Y[i]) == HUGE_VAL)
            throwError(d, "Coordinate Y at line %lu has caused overflow", keywordsLinesCount + i + 1);

        i++;
    }

    fclose(f);
    if (line) free(line);

    // print the coordinates data with enough log level
    for (size_t i = 0; i < d->nodesCount; i++)
        LOG(LOG_LVL_EVERYTHING, "Node %4lu: [%.2e, %.2e]", i+1, d->X[i], d->Y[i]);
}


static void getNameFromFile(char * line, int lineSize, char out[], Instance *d)
{
    memset(out, 0, 200); // set name array to 0

    char * nameBegin = strchr(line, ':');
    if (*(nameBegin+1) == ' ') 
        nameBegin += 2;

    int nameLen = lineSize + line - nameBegin - 1;
    if (nameLen > 200)
        throwError(d, "Name lenght exceeds 200 characters");
    memcpy(out, nameBegin, nameLen);
}

static void checkFileType(char * line, int lineSize, Instance *d)
{
    // here check if last 3 charaters(excludiing the '\n') are "TSP"
    char substr[3] = {0};
    memcpy(substr, &line[lineSize - 4], 3); // generate substring of 3 chars for logging/debugging purposes
    LOG(LOG_LVL_EVERYTHING, "Checking file TYPE keyword: comparing \"TSP\" with what is found at the of the line which is:%s", substr);
    if (strncmp(&line[lineSize - 4], "TSP", 3) != 0)
        throwError(d, "The file is either not of type TSP, or there are some characters (even blank spaces) after \"TSP\" and before the next line. \n\
                                     Check that the file used in input is of the correct type and correctly formatted. Only \"TSP\" files are currently supported");
}

static size_t getDimensionFromLine(char * line, int lineSize, Instance *d)
{
    // first find the pointer to the first number part of line and then convert it to integer checking for all errors
    char *numberFirstChar = strchr(line, ':');

    if (*(numberFirstChar + 1) == ' ')
        numberFirstChar += 2;
    else
        numberFirstChar++;
    // now numberFirstChar should be pointing to the first character that makes the decimal number in line
    LOG(LOG_LVL_EVERYTHING, "Getting the number of nodes from file: the first character of the number is:%c", *numberFirstChar);

    char *endPtr = NULL;
    size_t dimension = strtoul(numberFirstChar, &endPtr, 10);

    // check for errors on conversion
    if (endPtr == numberFirstChar)
        throwError(d, "Converting dimension number from file: first character that was supposed to be a number is not a number. \n\
                                        Dimension line in file is supposed to look like \"DIMENSION : <NUMBER>\" with just one separator \':\' and at most a \' \' after it before the numeric value");
    if (endPtr != &line[lineSize - 1])
        throwError(d, "Converting dimension number from file: there are unrecognized character before end of line with keyword DIMENSION");

    // check for error on converted number
    if (dimension == 0)
        throwError(d, "Could not properly convert dimension number in tsp file");
    if (dimension == ULONG_MAX)
        throwError(d, "Dimension value in file is too great. Either too much data(unlikely) or the dimension value is wrong");

    return dimension;
}

static size_t getEdgeWeightTypeFromLine(char * line, int lineSize, Instance *d)
{
    // first find the pointer to the weight type descriptor
    char *firstWgtTypePtr = strchr(line, ':');

    if (*(firstWgtTypePtr + 1) == ' ')
        firstWgtTypePtr += 2;
    else
        firstWgtTypePtr++;

    // now firstEdgeTypePtr should be pointing to the first character that describes the weight type
    LOG(LOG_LVL_EVERYTHING, "Getting the edge weight type from file: the first character of the weight type is:%c", *firstWgtTypePtr);

    int foundEdgeWeightTypeID = -1;
    for (size_t i = 0; i < EDGE_WEIGHT_TYPES_COUNT; i++)
    {
        size_t wgtTypeStrLen = strlen(wgtTypeStr[i]);

        if (strncmp(firstWgtTypePtr, wgtTypeStr[i], wgtTypeStrLen) == 0)
        {
            foundEdgeWeightTypeID = i;
            break;
        }
    }

    // check if no match has been found
    if (foundEdgeWeightTypeID == -1)
        throwError(d, "Getting the edge weight type from file: Could not identify the \"EDGE_WEIGTH_TYPE property\"");

    return foundEdgeWeightTypeID;
}

static inline int nProcessors()
{
    FILE *commandPipe;
    char *command = "nproc";
    char temp[10];
    commandPipe = (FILE*)popen(command, "r");
    fgets(temp, sizeof(temp), commandPipe);
    pclose(commandPipe);
    int numProcessors = atoi(temp);
    return numProcessors;
}

void saveSolution(Instance *d)
{  
    // create file for the solution
    char fileName[50] = "run/";
    strcat(fileName, d->params.name);
    strcat(fileName, ".opt.tour");
    FILE *solutionFile = fopen(fileName, "w");

    // inserting headers of the file
    fprintf(solutionFile, "NAME : %s\n", fileName);
    fprintf(solutionFile, "TYPE : TOUR\n");
    fprintf(solutionFile, "DIMENSION : %ld\n", d->nodesCount);
    fprintf(solutionFile, "TOUR_SECTION\n");

    // populating the file with the solution
    for (int i = 0; i < d->nodesCount; i++)
    {
        fprintf(solutionFile, "%d\n", d->solution.bestSolution[i] + 1);
    }
    fprintf(solutionFile, "-1\n");

    // closing the tour file
    fclose(solutionFile);
     
}

void plotSolution(Instance *d, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize)
{
    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"%s\"\n", d->params.name);
    fprintf(gnuplotPipe, "set terminal qt size %s\n", plotPixelSize);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pt 7 pointsize %d\n", pointColor, pointSize);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pointsize %d\n", tourPointColor, pointSize);

    // populating the plot
    
    fprintf(gnuplotPipe, "plot '-' with point linestyle 1, '-' with linespoint linestyle 2\n");

    // first plot only the points
    for (size_t i = 0; i < d->nodesCount; i++)
        fprintf(gnuplotPipe, "%f %f\n", d->X[i], d->Y[i]);
    fprintf(gnuplotPipe, "e\n");
    
    // second print the tour
    for (int i = 0; i <= d->nodesCount; i++)    // SOLUTION IS ALREDY SAVED IN ARRAY OF NODESCOUNT+1 ELEMENTS
    {
        int successorID = d->solution.bestSolution[i];
        fprintf(gnuplotPipe, "%f %f\n", d->X[successorID], d->Y[successorID]);
    }
    //fprintf(gnuplotPipe, "%lf %lf\n", d->X[d->solution.bestSolution[0]], d->Y[d->solution.bestSolution[0]]);
    fprintf(gnuplotPipe, "e\n");

    // force write on stream
    fflush(gnuplotPipe);

    // close stream
    pclose(gnuplotPipe);
}