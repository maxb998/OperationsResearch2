#ifndef _2OPT
#define _2OPT

#include "TspBase.h"

enum _2OptOptions
{
    _2OPT_BASE,
    _2OPT_AVX,
    _2OPT_PRECOMPUTED_COSTS
};

Solution _2OptBestFix(Solution *sol, enum _2OptOptions option);

double apply2OptBestFix(Solution *sol, enum _2OptOptions option);

#endif // _2OPT