#include "logLib.h"

static enum logLevel globLVL = LOG_LVL_DEBUG;

int LOG (enum logLevel lvl, char * line, ...)
{
    if (lvl < globLVL)
        return 0;

    va_list par;
    
    
    
}