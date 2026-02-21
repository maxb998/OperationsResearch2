#include "ArgParser.h"

#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


static void parserError (char * line, ...);

static void parseList(char *arg, const char separator, void *savePtr, int listLen, const char *paramName, void parserFunc());

static void parseUint(char *arg, char *expectedEndPtr, int *savePtr, int index, const char *paramName);

static void parseDouble(char *arg, char *expectedEndPtr, double *savePtr, int index, const char *paramName);

static void parseStringSubopts(char *arg, const char *subOptList[], int *savePtr, int subOptCount, const char *paramName);

static void printOptDoc(const char *doc, char *strBuff, int startPos, int terminalWidth);

static void printSuboptDoc(const char *doc, char *strBuff, int startPos, int terminalWidth);

static void printHelp(ArgOption *opt, ArgGroup *groups);


void parser(ArgOption *opt, ArgGroup *groups, int argc, char *argv[])
{
    // Data structure check
    for (int i = 0; groups[i].priority != -1; i++)
    {
        if (groups[i].priority < 0)
            parserError("groups must have all unique priorities greater than zero");

        for (int j = i+1; groups[j].priority != -1; j++)
            if (groups[i].priority == groups[j].priority)
                parserError("groups must have all unique priorities greater than zero");
    }

    // help msg
    for (int i = 1; i < argc; i++)
    {
        if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-?") == 0))
        {
            printHelp(opt, groups);
            exit(0);
        }
    }

    int optionsCount = 0;
    for (int i = 0; opt[i].key != -1; i++)
    {
        if (*(opt[i].name) == 0)
            parserError("option number %d does not have a name", i);
        if (opt[i].group == NULL)
            parserError("option number %d does not have a group", i);
        if (opt[i].count < 0)
            parserError("option \"%s\" has a listLen value of %d that is not valid", opt[i].name, opt[i].count);
        if (opt[i].subCount && (opt[i].dtype != DTYPE_STRING))
            parserError("option \"%s\" specifies a number of suboptions when its specified datatype is not string", opt[i].name);
        if ((opt[i].subCount < 0) || (opt[i].subCount == 1))
            parserError("option \"%s\" has a subCount value of %d that is not valid", opt[i].name, opt[i].subCount);
        if (opt[i].subCount && (opt[i].subNames == NULL))
            parserError("option \"%s\" has suboptions but no names were given", opt[i].name);
        if (opt[i].subCount && (opt[i].subDoc == NULL))
            parserError("option \"%s\" has suboptions but no documentation for them where provided", opt[i].name);
        if ((opt[i].dataPtr == NULL) && (opt[i].func == NULL))
            parserError("option \"%s\" has null dataPtr and func. This makes the option meaningless", opt[i].name);
        optionsCount++;
    }

    bool *requiredMap = malloc(optionsCount*2);
    if (requiredMap == NULL)
        parserError("Failed to allocate %d bytes of memory", optionsCount*2);
    bool *alreadyParsed = &(requiredMap[optionsCount]);
    for (int i = 0; opt[i].key != -1; i++)
    {
        requiredMap[i] = opt[i].required;
        alreadyParsed[i] = 0;
    }

    for (int i = 1; i < argc; i++) // start from 1 ignoring exec path
    {
        int optID = 0;
        if ((argv[i][0] == '-') && (argv[i][1] == '-'))
            while ((strcmp(opt[optID].name, &(argv[i][2])) != 0) && (opt[optID].key != -1))
                optID++;
        else if (argv[i][0] == '-')
            while ((opt[optID].key != argv[i][1]) && (opt[optID].key != -1))
                optID++;
        else // positional arg -> TODO
            parserError("positional arguments are not yet supported");

        requiredMap[optID] = false;
        if (alreadyParsed[optID])
            parserError("option \"--%s\"(or \"-%c\") can only be used once each execution", opt[optID].name, opt[optID].key);
        alreadyParsed[optID] = true;

        if (opt[optID].key == -1)
            parserError("option \"%s\" is not valid", argv[i]);

        int listLen = opt[optID].count;
        if (listLen == 0) listLen = 1;

        if (opt[optID].dataPtr)
        {
            switch (opt[optID].dtype)
            {
            case DTYPE_NONE:
                *(bool*)opt[optID].dataPtr = true;
                i--;
                break;
            case DTYPE_UINT:
                parseList(argv[i+1], ',', opt[optID].dataPtr, listLen, argv[i], parseUint);
                break;
            case DTYPE_DOUBLE:
                parseList(argv[i+1], ',', opt[optID].dataPtr, listLen, argv[i], parseDouble);
                break;
            case DTYPE_STRING:
                if (opt[optID].subCount == 0) // single string option
                {
                    char **ptr = opt[optID].dataPtr;
                    *ptr = argv[i+1];
                }
                else
                    parseStringSubopts(argv[i+1], opt[optID].subNames, opt[optID].dataPtr, opt[optID].subCount, argv[i]);
                break;
            }
        }

        i++;

        if (opt[optID].func)
            opt[optID].func(argv[i], opt[optID].dataPtr, opt[optID].funcData);
    }

    for (int i = 0; i < optionsCount; i++)
        if (requiredMap[i])
            parserError("missing required argument \"%s\"", opt[i].name);
    

    free(requiredMap);
}

static void parserError (char * line, ...)
{
    printf("ArgParser Error: ");

    va_list params;
    va_start(params, line);
    vprintf(line, params);
    va_end(params);

    printf("\n");

    exit(1);
}

static void parseList(char *arg, const char separator, void *savePtr, int listLen, const char *paramName, void parserFunc())
{
    char *endPtr = arg, *startPtr = arg;
    for (int i = 0; i < listLen; i++)
    {
        if (startPtr == NULL)
            parserError("missing and element for the option \"%s\". Check --help", paramName);
        
        while ((*endPtr != separator) && (*endPtr != 0))
            endPtr++;

        parserFunc(startPtr, endPtr, savePtr, i, paramName);

        endPtr++;
        startPtr = endPtr;
    }
}

static void parseUint(char *arg, char *expectedEndPtr, int *savePtr, int index, const char *paramName)
{
    char *endPtr;
    long num = strtol(arg, &endPtr, 10);
    if (num < 0)
        parserError("the value specified as \"%s\" cannot be negative", paramName);
    if (endPtr != expectedEndPtr)
        parserError("there are extra character after the \"%s\" value", paramName);

    savePtr[index] = num;
}

static void parseDouble(char *arg, char *expectedEndPtr, double *savePtr, int index, const char *paramName)
{
    char *endPtr;
    double num = strtod(arg, &endPtr);
    if (num <= 0)
        parserError("the value specified as \"%s\" must be a real number", paramName);
    if (endPtr != expectedEndPtr)
        parserError("there are extra character after the \"%s\" value", paramName);

    savePtr[index] = num;
}

static void parseStringSubopts(char *arg, const char *subOptList[], int *savePtr, int subOptCount, const char *paramName)
{
    for (int i = 0; i < subOptCount; i++)
    {
        if (strcmp(arg, subOptList[i]) == 0)
        {
            *savePtr = i;
            return;
        }
    }
    parserError("suboption \"%s\" specified with option \"%s\" was not recognized", arg, paramName);
}

static int printDocLineAndResetBuffer(char strBuff[], int writeIndex, int rstLen)
{
    strBuff[writeIndex] = 0;
    printf("%s\n", strBuff);
    
    writeIndex = 0;
    while (writeIndex < rstLen)
        strBuff[writeIndex++] = ' ';
    
    return writeIndex;
}

static void printOptDoc(const char *doc, char *strBuff, int startPos, int terminalWidth)
{
    int writeIndex = startPos;
    int i = 0;
    int doclen = strlen(doc);
    while (i < doclen)
    {
        // word detector
        int wordlen = 0;
        while ((i+wordlen < doclen) && (doc[i+wordlen] != ' '))
            wordlen++;
        
        // if word doesn't fit terminal width print and start new line
        if (writeIndex + wordlen > terminalWidth)
        {
            writeIndex = printDocLineAndResetBuffer(strBuff, writeIndex, startPos);
        }
        else
        {
            strncpy(&strBuff[writeIndex], &doc[i], wordlen+1);
            i += wordlen+1;
            writeIndex += wordlen+1;
        }
    }

    // print last word in doc string
    writeIndex = printDocLineAndResetBuffer(strBuff, writeIndex, terminalWidth);
}

static void printSuboptDoc(const char *doc, char *strBuff, int startPos, int terminalWidth)
{
    int writeIndex = startPos;
    int i = 0;
    int doclen = strlen(doc);
    while (i < doclen)
    {
        // word detector
        int wordlen = 0;
        while ((i+wordlen < doclen) && (doc[i+wordlen] != ' '))
            wordlen++;
        
        // if word doesn't fit terminal width print and start new line
        if (writeIndex + wordlen > terminalWidth)
        {
            writeIndex = printDocLineAndResetBuffer(strBuff, writeIndex, startPos);
        }
        else
        {
            strncpy(&strBuff[writeIndex], &doc[i], wordlen+1);
            i += wordlen+1;
            writeIndex += wordlen+1;
        }
    }

    // print last word in doc string
    writeIndex = printDocLineAndResetBuffer(strBuff, writeIndex, terminalWidth);
}

static void printHelp(ArgOption *opt, ArgGroup *groups)
{
    static const char *dtypeHelpStr[] = {
        " ",
        " <STRING> ",
        " <UINT> ",
        " <DOUBLE> "
    };

    // get terminal width
    unsigned short terminalWidth;
    {
        struct winsize w;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
        terminalWidth = w.ws_col;
    }

    int nGroups = 0;
    while (groups[nGroups].priority >= 0)
        nGroups++;

    ArgGroup **gsorted = malloc(sizeof(ArgGroup*)*nGroups);
    if (gsorted == NULL)
        parserError("failed to allocate memory for %u bytes", sizeof(ArgGroup*)*nGroups);

    for (int i = 0; i < nGroups; i++)
        gsorted[i] = &groups[i];

    // sort
    for (int i = 0; i < nGroups; i++)
    {
        for (int j = i+1; j < nGroups; j++)
        {
            if (gsorted[j]->priority < gsorted[i]->priority)
            {
                register ArgGroup *temp = gsorted[i];
                gsorted[i] = gsorted[j];
                gsorted[j] = temp;
            }
        }
    }

    char *strBuff = malloc(1024);
    if (strBuff == NULL)
        parserError("failed to allocate memory for %u bytes", 1024);

    for (int g = 0; g < nGroups; g++)
    {
        strncpy(strBuff, gsorted[g]->description, terminalWidth);
        printDocLineAndResetBuffer(strBuff, strlen(gsorted[g]->description), terminalWidth);

        // identify max length option to generate spacing
        int maxOptLen = 0;
        for (int i = 0; opt[i].key != -1; i++)
        {
            if (opt[i].group->priority != gsorted[g]->priority)
                continue;

            int currLen = 2;
            if (opt[i].key != 0)
                currLen += 3;
            currLen += strlen(opt[i].name) + 2;
            currLen += strlen(dtypeHelpStr[opt[i].dtype]);

            if (currLen > maxOptLen)
                maxOptLen = currLen;
        }
        if (maxOptLen == 0)
            parserError("internal error, maxOptLen == 0");

        for (int k = 0; opt[k].key != -1; k++)
        {
            if (opt[k].group->priority != gsorted[g]->priority)
                continue;

            int writeIndex = 2;
            if (opt[k].key != 0)
            {
                strBuff[writeIndex++] = '-';
                strBuff[writeIndex++] = opt[k].key;
                strBuff[writeIndex++] = ' ';
            }

            strBuff[writeIndex++] = '-';
            strBuff[writeIndex++] = '-';
            strcpy(&strBuff[writeIndex], opt[k].name);
            writeIndex += strlen(opt[k].name);
            
            strcpy(&strBuff[writeIndex], dtypeHelpStr[opt[k].dtype]);
            writeIndex += strlen(dtypeHelpStr[opt[k].dtype]);

            while (writeIndex < maxOptLen)
                strBuff[writeIndex++] = ' ';
            
            printOptDoc(opt[k].doc, strBuff, maxOptLen, terminalWidth);
            
            // print sub-options doc
            if (opt[k].subCount)
            {
                int maxSuboptLen = 0;
                for (int i = 0; i < opt[k].subCount; i++)
                {
                    int currLen = 3;
                    currLen += strlen(opt[k].subNames[i]) + 3;

                    if (currLen > maxSuboptLen)
                        maxSuboptLen = currLen;
                }

                for (int s = 0; s < opt[k].subCount; s++)
                {
                    writeIndex = maxOptLen;
                    strBuff[writeIndex++] = ' ';
                    strBuff[writeIndex++] = ' ';
                    strBuff[writeIndex++] = ' ';

                    strcpy(&strBuff[writeIndex], opt[k].subNames[s]);
                    writeIndex += strlen(opt[k].subNames[s]);
                    strBuff[writeIndex++] = ':';

                    while (writeIndex < maxOptLen + maxSuboptLen)
                        strBuff[writeIndex++] = ' ';

                    printSuboptDoc(opt[k].subDoc[s], strBuff, writeIndex, terminalWidth);
                }
                printf("\n");
            }
        }

        printf("\n");
    }

    free(strBuff);
    free(gsorted);
}
