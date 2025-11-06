#ifndef ARG_PARSER
#define ARGP_PARSER

enum ArgDType
{
    DTYPE_NONE, // flag type, no input
    DTYPE_STRING,
    DTYPE_UINT,
    DTYPE_DOUBLE
};

typedef struct
{
    const int priority; // priority of group over the other groups with -1 indicating the stop of group list. IT MUST BE UNIQUE TO EACH GROUP
    const char *description; // string that is printed at the beginnig of the group
} ArgGroup;

typedef struct
{
    const char key; // Key character that can be used to set the option at the execution of the program as "-KEY" (set to 0 if not used)
    const char *name; // Name of the option that is set when running the program and called with "--NAME"
    const enum ArgDType dtype; // Datatype of the option
    const ArgGroup *group; // pointer to group for grouping of options in help message
    void *dataPtr; // Pointer to the memory location in which save the parsed argument
    int count; // Number of arguments to parse (default:0 or 1 if dtype !NONE) (separator is ',')
    const char *doc; // Documentation string
    const char **subNames; // String array contaning the sub-options to match (only if dtype is string)
    const char **subDoc; // String array of documentation for each sub-option
    const unsigned int subCount; // Number of suboptions available
    void (*func)(); // Pointer to the function to call after standard parsing is done (function params must be (char *arg, void *dataPtr; void *funcData))
    void *funcData; // Pointer to data that will be given to the function called after the standard parsing
} ArgOption;

#define SIZE_CONST_ARR(arr) sizeof(arr)/sizeof(arr[0])

void parser(ArgOption *opt, ArgGroup *groups, int argc, char *argv[]);

#endif