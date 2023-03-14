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
    "TYPE",
    "DIMENSION",
    "EDGE_WEIGHT_TYPE",
    "NODE_COORD_SECTION"
};
#define KEYWORDS_COUNT 4

#define NO_KEYWORD_FOUND_ID -1
enum keywordID {
    TYPE_KEYWORD_ID, // 0
    DIMENSION_KEYWORD_ID, //1 
    EDGE_WEIGHT_TYPE_KEYWORD_ID, // 2
    NODE_COORD_SECTION_KEYWORD_ID // 3
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


int LOG (enum logLevel lvl, char * line, ...)
{
    // check log level
    if (lvl > globLVL) return 0;

    // get time
    /*time_t t;
    char dateStr[51];
    t = time(NULL);
    tzset();
    strftime(dateStr, sizeof(dateStr) - 1, "%a %b %d %T %Z %Y", localtime(&t));

    // print time
    printf("%s ", dateStr);*/

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
            {"avx", no_argument, 0, 'a'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int opt;

    while ((opt = getopt_long(argc, argv, "s:f:t:o:a", options, &option_index)) != -1)
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

        case 'a':
            d->params.useAVX = 1;
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
        LOG(LOG_LVL_DEBUG, line);

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
            LOG(LOG_LVL_ERROR, "Keyword \"%S\" is present more than one time. Check the .tsp file", keywords[keywordID]);

        switch (keywordID)
        {
        case TYPE_KEYWORD_ID:
            // here check if last 3 charaters(excludiing the '\n') are "TSP"
            char substr[3] = {0}; memcpy(substr, &line[lineSize - 4], 3);   // generate substring of 3 chars for logging/debugging purposes
            LOG(LOG_LVL_EVERYTHING, "Checking file TYPE keyword: comparing \"TSP\" with what is found at the of the line which is:%s", substr);
            if (strncmp(&line[lineSize - 4], "TSP", 3) != 0)
                LOG(LOG_LVL_ERROR, "The file is either not of type TSP, or there are some characters (even blank spaces) after \"TSP\" and before the next line. \n\
                                     Check that the file used in input is of the correct type and correctly formatted. Only \"TSP\" files are currently supported");
            // if this point is reached the has a correct type -> set the flag
            keywordsFound[TYPE_KEYWORD_ID] = 1;
            break;
        
        case DIMENSION_KEYWORD_ID:
            // first find the pointer to the first number part of line and then convert it to integer checking for all errors
            char * numberFirstChar = strchr(line, ':');

            if (*(numberFirstChar + 1) == ' ')
                numberFirstChar += 2;
            else
                numberFirstChar++;
            // now numberFirstChar should be pointing to the first character that makes the decimal number in line
            LOG(LOG_LVL_EVERYTHING, "Getting the number of nodes from file: the first character of the number is:%c", *numberFirstChar);

            char * endPtr = NULL;
            d->nodesCount = strtoul(numberFirstChar, &endPtr, 10);

            // check for errors on conversion
            if (endPtr == numberFirstChar)
                LOG(LOG_LVL_ERROR, "Converting dimension number from file: first character that was supposed to be a number is not a number. \n\
                                        Dimension line in file is supposed to look like \"DIMENSION : <NUMBER>\" with just one separator \':\' and at most a \' \' after it before the numeric value");
            if (endPtr != &line[lineSize-1])
                LOG(LOG_LVL_ERROR, "Converting dimension number from file: there are unrecognized character before end of line with keyword DIMENSION");

            // check for error on converted number
            if (d->nodesCount == 0)
                LOG(LOG_LVL_ERROR, "Could not properly convert dimension number in tsp file");
            if (d->nodesCount == ULONG_MAX)
                LOG(LOG_LVL_ERROR, "Dimension value in file is too great. Either too much data(unlikely) or the dimension value is wrong");

            // if this point is reached the file has a correct type -> set the flag
            keywordsFound[DIMENSION_KEYWORD_ID] = 1;
            break;

        case EDGE_WEIGHT_TYPE_KEYWORD_ID:
            // first find the pointer to the weight type descriptor
            char * firstWgtTypePtr = strchr(line, ':');

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
                LOG(LOG_LVL_ERROR, "Getting the edge weight type from file: Could not identify the \"EDGE_WEIGTH_TYPE property\"");
            
            d->params.edgeWeightType = foundEdgeWeightTypeID;

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
            LOG(LOG_LVL_ERROR, "Wierd error upon reading keywords from file");
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