#include "Tsp.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h> // needed to get the _POSIX_MONOTONIC_CLOCK and measure time

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
	"MAN_2D", // manhattan distance 2d
	"MAX_2D", // maximum distance 2d
    "EUC_2D", // euclidean distance 2d
	"CEIL_2D", // euclidean 2d rounded up
	"ATT", // special distance for problems att48 and att532
}; /*
	"EUC_3D", // euclidean distance 3d
	"MAN_3D", // manhattan distance 3d
	"MAX_3D", // maximum distance 3d
    "GEO", // geographical distance
	"XRAY1", // special distance for crystallography problems v1
	"XRAY2", // special distance for crystallography problems v2
	"EXPLICIT", // weights are specified in the file
	"SPECIAL" // special type of distance documented elsewhere
};// */
#define EDGE_WEIGHT_TYPES_COUNT 5

// file parsing functions

// gets the string specified in the "NAME" keyword in the file
static void getNameFromLine(char * line, int lineSize, char out[], Instance *inst);
// checks that the file type is "TSP"
static void checkFileType(char * line, int lineSize, Instance *inst);
// get the value related with the keyword "DIMENSION" and returns it as a int
static int getDimensionFromLine(char * line, int lineSize, Instance *inst);
// check that the string associated with "EDGE_WEIGHT_TYPE" is correct and return it as a number
static int getEdgeWeightTypeFromLine(char * line, int lineSize, Instance *inst);

static void readTspLine(char *line, int length, int index, Instance *inst, int keywordsLinesCount);

double readFile (Instance *inst)
{
    struct timespec timeStruct;
    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    double startTime = cvtTimespec2Double(timeStruct);

    FILE *f = fopen(inst->params.inputFile, "r");

    // check if was able to open file
    if (!f) { fprintf (stderr, "failed to open file for reading\n"); exit(EXIT_FAILURE); }

    char *line = NULL;
    size_t lineMemSize;
    int lineSize;
    int keywordsFound[KEYWORDS_COUNT] = { 0 };    // starting value assumes bad file

    int keywordsLinesCount = 0; // count the number of lines occupied by the keywords
    int keywordFinished = 0;
    while ((keywordFinished == 0) && ((lineSize = (int)getline(&line, &lineMemSize, f)) != EOF))
    {
        keywordsLinesCount++;
        //LOG(LOG_LVL_EVERYTHING, line);

        // check in the first part of each line if we find a useful keyword
        int keywordID = NO_KEYWORD_FOUND_ID;
        for (int i = 0; i < KEYWORDS_COUNT; i++)
        {
            int keywordLength = (int)strlen(keywords[i]);

            if (strncmp(line, keywords[i], keywordLength) == 0)
            {
                // if inside here means that keyword i has been found
                keywordID = i;
                break;
            }
        }

        if (keywordsFound[keywordID] == 1)
            throwError("Keyword \"%s\" is present more than one time. Check the .tsp file", keywords[keywordID]);

        switch (keywordID)
        {
        case NAME_KEYWORD_ID:
            getNameFromLine(line, lineSize, inst->params.name, inst);
            // set flag
            keywordsFound[NAME_KEYWORD_ID] = 1;
            break;

        case TYPE_KEYWORD_ID:
            checkFileType(line, lineSize, inst);
            // if this point is reached the has a correct type -> set the flag
            keywordsFound[TYPE_KEYWORD_ID] = 1;
            break;
        
        case DIMENSION_KEYWORD_ID:
            inst->nNodes = getDimensionFromLine(line, lineSize, inst);
            // set the flag
            keywordsFound[DIMENSION_KEYWORD_ID] = 1;
            break;

        case EDGE_WEIGHT_TYPE_KEYWORD_ID:
            inst->params.edgeWeightType  = getEdgeWeightTypeFromLine(line, lineSize, inst);
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
            throwError("Wierd error upon reading keywords from file. Check the code");
            break;
        }
    }

    // check all important keywords have been found
    for (int i = 0; i < KEYWORDS_COUNT; i++)
        if (keywordsFound[i] != 1)  // keyword has not been found in the loop above
            throwError("Important keyword \"%s\" has not been found/detected in the tsp file. Check the .tsp file", keywords[i]);

    // allocate memory
    // improve memory locality
    inst->X = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float)); // add AVX_VEC_SIZE extra elements at the end of each array for an easier time managing simd loads from this
    inst->Y = &inst->X[inst->nNodes + AVX_VEC_SIZE];

    // fill the memory with data
    int i = 0;
    while (((lineSize = getline(&line, &lineMemSize, f)) != EOF) && (strncmp(line, "EOF", 3) != 0))
    {
        // check if the number of data lines are not more than the number specified in dimension
        if (i >= inst->nNodes)
        {
            LOG(LOG_LVL_WARNING, "There may be more elements than the number specified with keyword \"DIMENSION\" in the .tsp file. Ignoring extra data");
            break;
        }

        readTspLine(line, lineSize, i, inst, keywordsLinesCount);

        i++;
    }

    if (i < inst->nNodes)
    {
        LOG(LOG_LVL_WARNING, "readFile: number of nodes in file are less than the number specified in the file. Setting total number of nodes to %ld", i);
        inst->nNodes = i;

        // reallocataing memory to better fit the smaller instance
        float *reallocX = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float));
        float *reallocY = &reallocX[inst->nNodes + AVX_VEC_SIZE];
        for (int j = 0; j < inst->nNodes; j++)
        {
            reallocX[j] = inst->X[j];
            reallocY[j] = inst->Y[j];
        }

        // free memory and swap pointers
        free(inst->X);
        inst->X = reallocX;
        inst->Y = reallocY;
    }

    // set to infinite the empty <AVX_VEC_SIZE> elements at the end of X and Y
    for (int i = inst->nNodes; i < inst->nNodes + AVX_VEC_SIZE; i++)
    {
        inst->X[i] = INFINITY;
        inst->Y[i] = INFINITY;
    }
    

    fclose(f);
    if (line) free(line);

    // print the coordinates data with enough log level
    //for (int i = 0; i < inst->nNodes; i++)
    //    LOG(LOG_LVL_EVERYTHING, "Node %4lu: [%.2e, %.2e]", i+1, inst->X[i], inst->Y[i]);

    clock_gettime(_POSIX_MONOTONIC_CLOCK, &timeStruct);
    return cvtTimespec2Double(timeStruct) - startTime;
}


static void getNameFromLine(char * line, int lineSize, char out[], Instance *inst)
{
    memset(out, 0, 200); // set name array to 0

    char * nameBegin = strchr(line, ':');
    if (*(nameBegin+1) == ' ') 
        nameBegin += 2;

    int nameLen = lineSize + line - nameBegin - 1;
    if (nameLen > 200)
        throwError("Name lenght exceeds 200 characters");
    memcpy(out, nameBegin, nameLen);
}

static void checkFileType(char * line, int lineSize, Instance *inst)
{
    // here check if last 3 charaters(excludiing the '\n') are "TSP"
    char substr[3] = {0};
    memcpy(substr, &line[lineSize - 4], 3); // generate substring of 3 chars for logging/debugging purposes
    //LOG(LOG_LVL_EVERYTHING, "Checking file TYPE keyword: comparing \"TSP\" with what is found at the of the line which is:%s", substr);
    if (strncmp(&line[lineSize - 4], "TSP", 3) != 0)
        throwError("The file is either not of type TSP, or there are some characters (even blank spaces) after \"TSP\" and before the next line. \nCheck that the file used in input is of the correct type and correctly formatted. Only \"TSP\" files are currently supported");
}

static int getDimensionFromLine(char * line, int lineSize, Instance *inst)
{
    // first find the pointer to the first number part of line and then convert it to integer checking for all errors
    char *numberFirstChar = strchr(line, ':');

    if (*(numberFirstChar + 1) == ' ')
        numberFirstChar += 2;
    else
        numberFirstChar++;
    // now numberFirstChar should be pointing to the first character that makes the decimal number in line
    //LOG(LOG_LVL_EVERYTHING, "Getting the number of nodes from file: the first character of the number is:%c", *numberFirstChar);

    char *endPtr = NULL;
    int dimension = (int)strtoul(numberFirstChar, &endPtr, 10);

    // check for errors on conversion
    if (endPtr == numberFirstChar)
        throwError("Converting dimension number from file: first character that was supposed to be a number is not a number. \n\
                                        Dimension line in file is supposed to look like \"DIMENSION : <NUMBER>\" with just one separator \':\' and at most a \' \' after it before the numeric value");
    if (endPtr != &line[lineSize - 1])
        throwError("Converting dimension number from file: there are unrecognized character before end of line with keyword DIMENSION");

    // check for error on converted number
    if (dimension == 0)
        throwError("Could not properly convert dimension number in tsp file");
    if (dimension == ULONG_MAX)
        throwError("Dimension value in file is too great. Either too much data(unlikely) or the dimension value is wrong");

    return dimension;
}

static int getEdgeWeightTypeFromLine(char * line, int lineSize, Instance *inst)
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
    for (int i = 0; i < EDGE_WEIGHT_TYPES_COUNT; i++)
    {
        int wgtTypeStrLen = (int)strlen(wgtTypeStr[i]);

        if (strncmp(firstWgtTypePtr, wgtTypeStr[i], wgtTypeStrLen) == 0)
        {
            foundEdgeWeightTypeID = i;
            break;
        }
    }

    // check if no match has been found
    if (foundEdgeWeightTypeID == -1)
        throwError("Getting the edge weight type from file: Could not identify the \"EDGE_WEIGTH_TYPE property\"");

    return foundEdgeWeightTypeID;
}

static void readTspLine(char *line, int length, int index, Instance *inst, int keywordsLinesCount)
{
    double numbers[3] = { 0 };
    char *endPtr = line;

    int i = 0, numID = 0;
    while (i < length)
    {
        // scroll past any black spaces
        if ((line[i] == ' ') || (line[i] == '\n'))
        {
            i++;
            continue;
        }

        if (numID >= 3)
            throwError("Line %lu is not formatted correctly. More than 3 numbers for each line are not supported", keywordsLinesCount + index + 1);

        // convert number and save into array
        char *oldEndPtr = endPtr;
        numbers[numID] = strtod(&line[i], &endPtr);
        if (endPtr == oldEndPtr) // error on strtod
            throwError("Conversion at line %lu has gone wrong", keywordsLinesCount + index + 1);
        numID++;

        // move i to end of number
        i += endPtr - oldEndPtr;
    }

    if ((int)numbers[0] != index + 1)
        LOG(LOG_LVL_WARNING, "readFile: File has inconsistent enumeration, node %lu at line %lu should be identified with %ld", (int)numbers[0], keywordsLinesCount + index + 1, index + 1);

    inst->X[index] = numbers[1];
    inst->Y[index] = numbers[2];
}

// cmd flags that are not meaningful to the output of the run, therefore not worth printing in the file with the solution
#define MEANINGLESS_FLAGS_COUNT 4
const char *meaninglessFlags[] = {
    "--plot",
    "-p",
    "--save",
    "-s"
};
// cmd params that are not meaningful to the output of the run, therefore not worth printing in the file with the solution
#define MEANINGLESS_ARGS_COUNT 4
const char *meaninglessArgs[] = {
    "--loglvl",
    "-l",
    "--threads",
    "-j"
};

void saveSolution(Solution *sol, int argc, char *argv[])
{
    // define file name
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char datetime[40] = { 0 }; // 21 should be enough but gives warning in compilation and I don't like warnings
    sprintf(datetime, "_%d-%02d-%02d_%02d:%02d:%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    
    char fileName[300] = { 0 };
    strcat(fileName, "run/");
    strcat(fileName, sol->instance->params.name);
    strcat(fileName, datetime);
    strcat(fileName, ".tour");

    FILE *solutionFile = fopen(fileName, "w");

    // inserting headers of the file
    fprintf(solutionFile, "NAME : %s.tour\n", sol->instance->params.name);

    // add description of the parameters used for this run
    char comment[1000] = { 0 };
    for (int i = 1; i < argc-1; i++)
    {
        // check if args is meaningless to save
        int skip = 0;
        for (int j = 0; j < MEANINGLESS_FLAGS_COUNT; j++)
            if (strcmp(argv[i], meaninglessFlags[j]) == 0)
                { skip = 1; break; }
        for (int j = 0; j < MEANINGLESS_ARGS_COUNT; j++)
            if (strcmp(argv[i], meaninglessArgs[j]) == 0)
                { i++; skip = 1; break; }
        if (skip) continue;

        // avoid saving the dashes
        char *argStart = argv[i];
        for (; *argStart == '-'; argStart++){}
        
        // save the args
        if ((i + 1 < argc) && (argv[i + 1][0] != '-'))
        {
            strcat(comment, argStart);
            strcat(comment, "=");
            i++;
            strcat(comment, argv[i]);
        }
        else
            strcat(comment, argStart);
        
        strcat(comment, "  ");
    }
    fprintf(solutionFile, "COMMENT : %s\n", comment);

    fprintf(solutionFile, "TYPE : TOUR\n");
    fprintf(solutionFile, "DIMENSION : %d\n", sol->instance->nNodes);
    fprintf(solutionFile, "TOUR_SECTION\n");

    // populating the file with the solution
    for (int i = 0; i < sol->instance->nNodes; i++)
        fprintf(solutionFile, "%d\n", sol->indexPath[i] + 1);
    fprintf(solutionFile, "-1\n");

    // closing the tour file
    fclose(solutionFile);
}

void plotSolution(Solution *sol, const char * plotPixelSize, const char * pointColor, const char * tourPointColor, const int pointSize, const bool printIndex)
{
    Instance *inst = sol->instance;
    int n = inst->nNodes;

    // creating the pipeline for gnuplot
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    // gnuplot settings
    fprintf(gnuplotPipe, "set title \"%s\"\n", sol->instance->params.name);
    fprintf(gnuplotPipe, "set terminal qt size %s\n", plotPixelSize);

    // set plot linestyles
    fprintf(gnuplotPipe, "set style line 1 linecolor rgb '%s' pt 7 pointsize %d\n", pointColor, pointSize);
    fprintf(gnuplotPipe, "set style line 2 linecolor rgb '%s' pointsize 0\n", tourPointColor);//, pointSize);


    // assign number to points
    if (printIndex)
        for (int i = 0; i < n; i++)
            fprintf(gnuplotPipe, "set label \"%d\" at %f,%f\n", i, inst->X[i], inst->Y[i]);

    // populating the plot
    
    fprintf(gnuplotPipe, "plot '-' with point linestyle 1, '-' with linespoint linestyle 2\n");

    // first plot only the points
    for (int i = 0; i < n; i++)
        fprintf(gnuplotPipe, "%f %f\n", inst->X[i], inst->Y[i]);
    fprintf(gnuplotPipe, "e\n");

    // second print the tour
    for (int i = 0; i < n; i++)
        fprintf(gnuplotPipe, "%f %f\n", inst->X[sol->indexPath[i]], inst->Y[sol->indexPath[i]]);
    fprintf(gnuplotPipe, "%f %f\n", inst->X[sol->indexPath[0]], inst->Y[sol->indexPath[0]]);
    
    fprintf(gnuplotPipe, "e\n");

    // force write on stream
    fflush(gnuplotPipe);

    // close stream
    pclose(gnuplotPipe);
}