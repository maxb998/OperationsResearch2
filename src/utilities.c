#include "utilities.h"

static enum logLevel globLVL = LOG_LVL_DEBUG;

const char * logLevelString [] = {
	"NONE", // 0
	"CRIT", // 1
	"WARN", // 2
	"NOTI", // 3
	"LOG ", // 4
	"DEBG", // 5
	"ALL"	// 6
};

int LOG (enum logLevel lvl, char * line, ...)
{
    // check log level
    if (lvl > globLVL) return 0;

    // get time
    time_t t;
    char dateStr[51];
    t = time(NULL);
    tzset();
    strftime(dateStr, sizeof(dateStr) - 1, "%a %b %d %T %Z %Y", localtime(&t));

    // print level of log and time
    printf("%s [%s] : ", dateStr, logLevelString[lvl]);

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
            {
                printf("ERROR: File \"%s\" not found\n", optarg);
                exit(EXIT_FAILURE);
            }
            strncpy(d->params.inputFile, optarg, strlen(optarg)+1);
            break;
        
        case 't':
            d->params.threadsCount = atoi(optarg);
            break;

        case 'a':
            d->params.useAVX = 1;
            break;
        
        default:
            abort();
        }
    }
    
    LOG(LOG_LVL_DEBUG, "Received arguments:");
    LOG(LOG_LVL_DEBUG,"    Random Seed  = %d", d->params.randomSeed);
    LOG(LOG_LVL_DEBUG,"    Filename     = %s", d->params.inputFile);
    LOG(LOG_LVL_DEBUG,"    Thread Count = %d", d->params.threadsCount);
}

void readFile (Instance *d)
{
    FILE *f = fopen(d->params.inputFile, "r");

    // check if was able to open file
    if (!f) { fprintf (stderr, "failed to open file for reading\n"); exit(EXIT_FAILURE); }

    char *line = NULL;
    size_t lineMemSize;
    int lineSize;
    int keywordFound[] = { 0, 0, 0, 0 };    // starting value assumes bad file

    while ((lineSize = getline(&line, &lineMemSize, f)) != EOF)
    {
        LOG(LOG_LVL_EVERYTHING, line);

        if (strstr(line, "TYPE"))   // check type
        {
            if (!strstr(line, "TSP") && !strstr(line, "tsp"))   // check for both upper and lower case NO MIXED TYPES 
            {
                LOG(LOG_LVL_CRITICAL, "ERROR: File type is not TSP");
                exit(EXIT_FAILURE);
            }
            if (keywordFound[0] == 1)
            {
                LOG(LOG_LVL_CRITICAL, "ERROR: File contains more than one \"TYPE\" keyword");
                exit(EXIT_FAILURE);
            }
            keywordFound[0] = 1;  // remember that first keyword has been found and file type is correct
        }
        else if (strstr(line, "DIMENSION")) // get dimension
        {
            // get pointer to beginning of dimesion number
            char * substrBeginPtr = strchr(line, ':');
            if (!substrBeginPtr)    // exit with error if no separator ':' is found
            {
                LOG(LOG_LVL_CRITICAL, "ERROR: File \"DIMENSION\" line is missing separator \':\'");
                exit(EXIT_FAILURE);
            }

            // check if there is a ' ' after ':' and, if there is a ' ' move the pointer forward by 2 (need to point to the first char of the number)
            if (*(substrBeginPtr+1) == ' ') substrBeginPtr += 2;

            // get length of substring containing only the dimension number
            int substrLen = line + strlen(line) - substrBeginPtr - 1;   // assumes that lines ends with '\n'

            // create substring with the value of "DIMENSION" only
            char substr[20] = { 0 };
            memcpy(substr, substrBeginPtr, substrLen);

            LOG(LOG_LVL_EVERYTHING, "TSP file reading -> Get Dimension -> substr = \"%s\"", substr);

            // convert substring to nodeCount and check for errors
            char * endPtr = NULL;
            d->nodesCount = strtod(substr, &endPtr);
            if (endPtr != &substr[substrLen])   // check if all substring was converted. if not throw error and exit
            {
                LOG(LOG_LVL_CRITICAL, "ERROR: \"DIMENSION\" keyword value conversion to double has gone wrong");
                exit(EXIT_FAILURE);
            }

        }
        
    }
    exit(0);

    // allocate memory
    d->coords = malloc(d->nodesCount * 2 * sizeof(double));

    // fill the matrix
    int placeholder = 0;
    for (size_t i = 0; i < d->nodesCount; i++)
    {
        fscanf(f, "%d %lf %lf", &placeholder, &d->coords[i], &d->coords[i + d->nodesCount]);

        // print acquired data
        LOG(LOG_LVL_EVERYTHING, "Node %d : [%6.3lf, %6.3lf]", placeholder, d->coords[i], d->coords[i + d->nodesCount]);
    }
    
    
    if(line) free(line);

    fclose(f);
}