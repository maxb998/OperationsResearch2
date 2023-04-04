#ifndef _2OPT
#define _2OPT

#include "TspBase.h"

enum _2OptOptions
{
    _2OPT_AVX_ST,
    _2OPT_BASE_MT,
    _2OPT_AVX_MT,
    _2OPT_PRECOMPUTED_COSTS_MT
};

Solution _2OptBestFix(Solution *sol, enum _2OptOptions option);

double apply2OptBestFix(Solution *sol, enum _2OptOptions option);

#endif // _2OPT