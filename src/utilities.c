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
	"EUC_2D", // 0, euclidean distance 2d
	"MAN_2D", // 1, manhattan distance 2d
	"MAX_2D", // 2, maximum distance 2d
	"CEIL_2D", // 3, euclidean 2d rounded up
    "GEO", // 4,
	"ATT", // 5, special distance for problems att48 and att532
	"XRAY1", // 6, special distance for crystallography problems v1
	"XRAY2", // 7, special distance for crystallography problems v2
	"EUC_3D", // 8, euclidean distance 3d
	"MAN_3D", // 9, manhattan distance 3d
	"MAX_3D", // 10, maximum distance 3d
	"EXPLICIT", // 11, weights are specified in the file
	"SPECIAL" // 12, special type of distance documented elsewhere
};
#define EDGE_WEIGHT_TYPES_COUNT 12

// file parsing functions
void getNameFromFile(char * line, int lineSize, char out[]);
void checkFileType(char * line, int lineSize);
size_t getDimensionFromLine(char * line, int lineSize);
size_t getEdgeWeightTypeFromLine(char * line, int lineSize);


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

    if (lvl == LOG_LVL_ERROR)
        exit(EXIT_FAILURE);

    return 0;
}

void parseArgs (Instance *d, int argc, char *argv[])
{
    static struct option options[] = 
        {
            {"seed", required_argument, 0, 's'},
            {"s", required_argument, 0, 's'},
            {"file", required_argument, 0, 'f'},
            {"f", required_argument, 0, 'f'},
            {"threads", required_argument, 0, 't'},
            {"t", required_argument, 0, 't'},
            {"out", required_argument, 0, 'o'},
            {"o", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int opt;

    while ((opt = getopt_long(argc, argv, "s:f:t:o:", options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 's':
            d->params.randomSeed = atof(optarg);
            break;

        case 'f':
            if (access(optarg, R_OK) != 0)
                LOG(LOG_LVL_ERROR, "ERROR: File \"%s\" not found\n", optarg);

            strncpy(d->params.inputFile, optarg, strlen(optarg)+1);
            break;
        
        case 't':
            d->params.threadsCount = strtoul(optarg, NULL, 10);
            break;
        
        default:
            abort();
        }
    }
    
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
            LOG(LOG_LVL_ERROR, "Keyword \"%s\" is present more than one time. Check the .tsp file", keywords[keywordID]);

        switch (keywordID)
        {
        case NAME_KEYWORD_ID:
            getNameFromFile(line, lineSize, d->params.name);
            // set flag
            keywordsFound[NAME_KEYWORD_ID] = 1;
            break;

        case TYPE_KEYWORD_ID:
            checkFileType(line, lineSize);
            // if this point is reached the has a correct type -> set the flag
            keywordsFound[TYPE_KEYWORD_ID] = 1;
            break;
        
        case DIMENSION_KEYWORD_ID:
            d->nodesCount = getDimensionFromLine(line, lineSize);
            // set the flag
            keywordsFound[DIMENSION_KEYWORD_ID] = 1;
            break;

        case EDGE_WEIGHT_TYPE_KEYWORD_ID:
            d->params.edgeWeightType = getEdgeWeightTypeFromLine(line, lineSize);
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
            LOG(LOG_LVL_ERROR, "Wierd error upon reading keywords from file. Check the code");
            break;
        }
    }

    // check all important keywords have been found
    for (size_t i = 0; i < KEYWORDS_COUNT; i++)
        if (keywordsFound[i] != 1)  // keyword has not been found in the loop above
            LOG(LOG_LVL_ERROR, "Important keyword \"%s\" has not been found/detected in the tsp file. Check the .tsp file", keywords[i]);

    // allocate memory
    size_t memElemsToAlloc = d->nodesCount + AVX_VEC_SIZE - d->nodesCount % AVX_VEC_SIZE;
    d->X = aligned_alloc(32, memElemsToAlloc * sizeof(double));
    d->Y = aligned_alloc(32, memElemsToAlloc * sizeof(double));

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
            LOG(LOG_LVL_ERROR, "Conversion at line %lu of the line number has gone wrong(must not be 0)", keywordsLinesCount + i + 1);
        
        char * separatorPtr = strchr(line, ' ');
        if (!separatorPtr)
            LOG(LOG_LVL_ERROR, "Space separator at line %lu of file has not been found");
        if (endPtr != separatorPtr)
            LOG(LOG_LVL_ERROR, "Conversion at line %lu of the line number has gone wrong(must be an integer separated from other value by a space at the end)", keywordsLinesCount + i + 1);
        
        if (lineNumber != i+1)
            LOG(LOG_LVL_ERROR, "Line number at line %lu of the file is not what is supposed to be (%lu instead of %lu)", keywordsLinesCount + i + 1, lineNumber, i+1);

        // now we can convert the actual coordinates
        char * xCoordStrPtr = separatorPtr + 1;
        separatorPtr = strchr(xCoordStrPtr, ' ');
        if (!separatorPtr)
            LOG(LOG_LVL_ERROR, "Space separator at line %lu of file has not been found");

        d->X[i] = strtod(xCoordStrPtr, &endPtr);

        if (xCoordStrPtr == endPtr)
            LOG(LOG_LVL_ERROR, "Conversion of X coordinate at line %lu has gone wrong. Check the .tsp file", keywordsLinesCount + i + 1);
        if (endPtr != separatorPtr)
            LOG(LOG_LVL_ERROR, "Conversion of X coordinate at line %lu has gone wrong. There are unwanted characters at the end ", keywordsLinesCount + i + 1);
        if (d->X[i] == HUGE_VAL)
            LOG(LOG_LVL_ERROR, "Coordinate X at line %lu has caused overflow", keywordsLinesCount + i + 1);
        if (d->X[i] == __DBL_MIN__)
            LOG(LOG_LVL_ERROR, "Coordinate X at line %lu has caused underflow", keywordsLinesCount + i + 1);

        char * yCoordStrPtr = separatorPtr + 1;
        separatorPtr = strchr(yCoordStrPtr, '\n');
        if (!separatorPtr)
            LOG(LOG_LVL_ERROR, "'\n' at line %lu of file has not been found");
        
        d->Y[i] = strtod(yCoordStrPtr, &endPtr);

        if (yCoordStrPtr == endPtr)
            LOG(LOG_LVL_ERROR, "Conversion of Y coordinate at line %lu has gone wrong. Check the .tsp file", keywordsLinesCount + i + 1);
        if (endPtr != separatorPtr)
            LOG(LOG_LVL_ERROR, "Conversion of Y coordinate at line %lu has gone wrong. There are unwanted characters at the end ", keywordsLinesCount + i + 1);
        if (d->Y[i] == HUGE_VAL)
            LOG(LOG_LVL_ERROR, "Coordinate Y at line %lu has caused overflow", keywordsLinesCount + i + 1);
        if (d->Y[i] == __DBL_MIN__)
            LOG(LOG_LVL_ERROR, "Coordinate Y at line %lu has caused underflow", keywordsLinesCount + i + 1);

        i++;
    }

    fclose(f);
    if (line) free(line);

    // print the coordinates data with enough log level
    for (size_t i = 0; i < d->nodesCount; i++)
        LOG(LOG_LVL_DEBUG, "Node %3lu: [%.2lf, %.2lf]", i+1, d->X[i], d->Y[i]);
}

void saveSolution(Instance *d)
{  
    // create file for the solution
    char fileName[50] = "run/demo.opt.tour";
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

void plotSolution(Instance *d)
{
    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    char fileName[50] = "DEMO PLOT";            // DA INSERIRE IL VERO NOME DEL FILE
    fprintf(gnuplotPipe, "set title \"%s\"\n", fileName);

    // populating the plot
    fprintf(gnuplotPipe, "plot '-' \n");
    int i = 0;
    do
    {
      fprintf(gnuplotPipe, "%lf %lf\n", d->X[i], d->Y[i]);
      printf("%lf %lf\n", d->X[i], d->Y[i]);
      i = d->solution.bestSolution[i];
    } while (i != 0);
    fprintf(gnuplotPipe, "e\n");
    fflush(gnuplotPipe);
}