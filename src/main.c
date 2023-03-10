#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>

#include "utilities.h"
#include "tsp.h"

int main (int argc, char *argv[])
{
    Instance d; d.params.useAVX = 0;    // avx disabled by default
    parseArgs(&d, argc, argv);
    readFile(&d);

    return EXIT_SUCCESS;
}