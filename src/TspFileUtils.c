#include "TspFileUtils.h"
#include "TspUtilities.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <math.h>

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
static void getNameFromFile(char * line, int lineSize, char out[], Instance *inst);
// checks that the file type is "TSP"
static void checkFileType(char * line, int lineSize, Instance *inst);
// get the value related with the keyword "DIMENSION" and returns it as a size_t
static size_t getDimensionFromLine(char * line, int lineSize, Instance *inst);
// check that the string associated with "EDGE_WEIGHT_TYPE" is correct and return it as a number
static size_t getEdgeWeightTypeFromLine(char * line, int lineSize, Instance *inst);


void readFile (Instance *inst)
{
    FILE *f = fopen(inst->params.inputFile, "r");

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
            throwError(inst, NULL, "Keyword \"%s\" is present more than one time. Check the .tsp file", keywords[keywordID]);

        switch (keywordID)
        {
        case NAME_KEYWORD_ID:
            getNameFromFile(line, lineSize, inst->params.name, inst);
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
            inst->params.edgeWeightType = getEdgeWeightTypeFromLine(line, lineSize, inst);
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
            throwError(inst, NULL, "Wierd error upon reading keywords from file. Check the code");
            break;
        }
    }

    // check all important keywords have been found
    for (size_t i = 0; i < KEYWORDS_COUNT; i++)
        if (keywordsFound[i] != 1)  // keyword has not been found in the loop above
            throwError(inst, NULL, "Important keyword \"%s\" has not been found/detected in the tsp file. Check the .tsp file", keywords[i]);

    // allocate memory
    // improve memory locality
    inst->X = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float)); // add AVX_VEC_SIZE extra elements at the end of each array for an easier time managing simd loads from this
    inst->Y = &inst->X[inst->nNodes + AVX_VEC_SIZE];

    // fill the memory with data
    size_t i = 0;
    while (((lineSize = getline(&line, &lineMemSize, f)) != EOF) && (strncmp(line, "EOF", 3) != 0))
    {
        // check if the number of data lines are not more than the number specified in dimension
        if (i >= inst->nNodes)
        {
            LOG(LOG_LVL_WARNING, "There may be more elements than the number specified with keyword \"DIMENSION\" in the .tsp file. Ignoring extra data");
            break;
        }

        // get the line number, convert it and check it's ok
        char * endPtr = NULL;
        size_t lineNumber = strtoul(line, &endPtr, 10);
        if (lineNumber == 0)
            throwError(inst, NULL, "Conversion at line %lu of the line number has gone wrong(must not be 0)", keywordsLinesCount + i + 1);
        else if (lineNumber != i+1)
            LOG(LOG_LVL_WARNING, "readFile: File has inconsistent enumeration, node %ld at line %ld should be identified with %ld", lineNumber, keywordsLinesCount + i + 1, i+1);
        
        
        char * separatorPtr = strchr(line, ' ');
        if (!separatorPtr)
            throwError(inst, NULL, "Space separator at line %lu of file has not been found");
        if (endPtr != separatorPtr)
            throwError(inst, NULL, "Conversion at line %lu of the line number has gone wrong(must be an integer separated from other value by a space at the end)", keywordsLinesCount + i + 1);

        // now we can convert the actual coordinates
        char * xCoordStrPtr = separatorPtr + 1;
        separatorPtr = strchr(xCoordStrPtr, ' ');
        if (!separatorPtr)
            throwError(inst, NULL, "Space separator at line %lu of file has not been found");

        inst->X[i] = strtof(xCoordStrPtr, &endPtr);

        if (xCoordStrPtr == endPtr)
            throwError(inst, NULL, "Conversion of X coordinate at line %lu has gone wrong. Check the .tsp file", keywordsLinesCount + i + 1);
        if (endPtr != separatorPtr)
            throwError(inst, NULL, "Conversion of X coordinate at line %lu has gone wrong. There are unwanted characters at the end ", keywordsLinesCount + i + 1);
        if (fabsf(inst->X[i]) == HUGE_VALF)
            throwError(inst, NULL, "Coordinate X at line %lu has caused overflow", keywordsLinesCount + i + 1);

        char * yCoordStrPtr = separatorPtr + 1;
        separatorPtr = strchr(yCoordStrPtr, '\n');
        if (!separatorPtr)
            throwError(inst, NULL, "'\n' at line %lu of file has not been found");
        
        inst->Y[i] = strtof(yCoordStrPtr, &endPtr);

        if (yCoordStrPtr == endPtr)
            throwError(inst, NULL, "Conversion of Y coordinate at line %lu has gone wrong. Check the .tsp file", keywordsLinesCount + i + 1);
        if (endPtr != separatorPtr)
            throwError(inst, NULL, "Conversion of Y coordinate at line %lu has gone wrong. There are unwanted characters at the end ", keywordsLinesCount + i + 1);
        if (fabsf(inst->Y[i]) == HUGE_VAL)
            throwError(inst, NULL, "Coordinate Y at line %lu has caused overflow", keywordsLinesCount + i + 1);

        i++;
    }

    if (i < inst->nNodes)
    {
        LOG(LOG_LVL_WARNING, "readFile: number of nodes in file are less than the number specified in the file. Setting total number of nodes to %ld", i);
        inst->nNodes = i;

        // reallocataing memory to better fit the smaller instance
        float *reallocX = malloc((inst->nNodes + AVX_VEC_SIZE) * 2 * sizeof(float));
        float *reallocY = &reallocX[inst->nNodes + AVX_VEC_SIZE];
        for (size_t j = 0; j < inst->nNodes; j++)
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
    for (size_t i = inst->nNodes; i < inst->nNodes + AVX_VEC_SIZE; i++)
    {
        inst->X[i] = INFINITY;
        inst->Y[i] = INFINITY;
    }
    

    fclose(f);
    if (line) free(line);

    // print the coordinates data with enough log level
    for (size_t i = 0; i < inst->nNodes; i++)
        LOG(LOG_LVL_EVERYTHING, "Node %4lu: [%.2e, %.2e]", i+1, inst->X[i], inst->Y[i]);
}


static void getNameFromFile(char * line, int lineSize, char out[], Instance *inst)
{
    memset(out, 0, 200); // set name array to 0

    char * nameBegin = strchr(line, ':');
    if (*(nameBegin+1) == ' ') 
        nameBegin += 2;

    int nameLen = lineSize + line - nameBegin - 1;
    if (nameLen > 200)
        throwError(inst, NULL, "Name lenght exceeds 200 characters");
    memcpy(out, nameBegin, nameLen);
}

static void checkFileType(char * line, int lineSize, Instance *inst)
{
    // here check if last 3 charaters(excludiing the '\n') are "TSP"
    char substr[3] = {0};
    memcpy(substr, &line[lineSize - 4], 3); // generate substring of 3 chars for logging/debugging purposes
    LOG(LOG_LVL_EVERYTHING, "Checking file TYPE keyword: comparing \"TSP\" with what is found at the of the line which is:%s", substr);
    if (strncmp(&line[lineSize - 4], "TSP", 3) != 0)
        throwError(inst, NULL, "The file is either not of type TSP, or there are some characters (even blank spaces) after \"TSP\" and before the next line. \n\
                                     Check that the file used in input is of the correct type and correctly formatted. Only \"TSP\" files are currently supported");
}

static size_t getDimensionFromLine(char * line, int lineSize, Instance *inst)
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
        throwError(inst, NULL, "Converting dimension number from file: first character that was supposed to be a number is not a number. \n\
                                        Dimension line in file is supposed to look like \"DIMENSION : <NUMBER>\" with just one separator \':\' and at most a \' \' after it before the numeric value");
    if (endPtr != &line[lineSize - 1])
        throwError(inst, NULL, "Converting dimension number from file: there are unrecognized character before end of line with keyword DIMENSION");

    // check for error on converted number
    if (dimension == 0)
        throwError(inst, NULL, "Could not properly convert dimension number in tsp file");
    if (dimension == ULONG_MAX)
        throwError(inst, NULL, "Dimension value in file is too great. Either too much data(unlikely) or the dimension value is wrong");

    return dimension;
}

static size_t getEdgeWeightTypeFromLine(char * line, int lineSize, Instance *inst)
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
        throwError(inst, NULL, "Getting the edge weight type from file: Could not identify the \"EDGE_WEIGTH_TYPE property\"");

    return foundEdgeWeightTypeID;
}

void saveSolution(Solution *sol)
{  
    // create file for the solution
    char fileName[50] = "run/";
    strcat(fileName, sol->instance->params.name);
    strcat(fileName, ".opt.tour");
    FILE *solutionFile = fopen(fileName, "w");

    // inserting headers of the file
    fprintf(solutionFile, "NAME : %s\n", fileName);
    fprintf(solutionFile, "TYPE : TOUR\n");
    fprintf(solutionFile, "DIMENSION : %ld\n", sol->instance->nNodes);
    fprintf(solutionFile, "TOUR_SECTION\n");

    int *solID = sol->indexPath;

    // populating the file with the solution
    for (int i = 0; i < sol->instance->nNodes; i++)
        fprintf(solutionFile, "%d\n", solID[i] + 1);
    fprintf(solutionFile, "-1\n");

    

    // closing the tour file
    fclose(solutionFile);
    free(solID);
}