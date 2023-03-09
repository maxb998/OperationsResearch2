#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>

#include "tsp.h"

void parseArgs (Instance *d, int argc, char *argv[]);
void readFile (Instance *d);

int main (int argc, char *argv[])
{
    Instance d; d.params.useAVX = 0;    // avx disabled by default
    parseArgs(&d, argc, argv);
    readFile(&d);

    return EXIT_SUCCESS;
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
            strncpy(d->params.inputFile, optarg, strlen(optarg));
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
    
    if (VERBOSE >= 70)
    {
        printf("Received arguments:\n");
        printf("    Random Seed  = %d\n", d->params.randomSeed);
        printf("    Filename     = %s\n", d->params.inputFile);
        printf("    Thread Count = %d\n", d->params.threadsCount);
        printf("    AVX Enabled  = ");  if (d->params.useAVX == 0) printf("NO\n"); else printf("YES\n");
        printf("\n");
    }
}

void readFile (Instance *d)
{
    FILE *f = fopen(d->params.inputFile, "r");

    // check if was able to open file
    if (!f) { fprintf (stderr, "failed to open file for reading\n"); exit(EXIT_FAILURE); }

    char *line = NULL;
    size_t lineMemSize;
    int lineSize;

    for (size_t i = 0; (i < 6) && ((lineSize = getline(&line, &lineMemSize, f)) != EOF); i++)
    {
        if (VERBOSE >= 100) printf("%s", line); // verobse option

        if (i == 3)
        {
            sscanf(line, "DIMENSION : %lu", &d->nodesCount);
            //printf("%lu\n", d->nodesCount);
        }
    }

    // allocate memory
    d->coords = malloc(d->nodesCount * 2 * sizeof(double));

    // fill the matrix
    int placeholder = 0;
    for (size_t i = 0; i < d->nodesCount; i++)
    {
        fscanf(f, "%d %lf %lf", &placeholder, &d->coords[i], &d->coords[i + d->nodesCount]);

        // print acquired data
        if (VERBOSE >= 150) printf("Node %d : [%6.3lf, %6.3lf]\n", placeholder, d->coords[i], d->coords[i + d->nodesCount]);
    }
    
    
    if(line) free(line);

    fclose(f);
}