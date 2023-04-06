#ifndef EXTRA_MILEAGE
#define EXTRA_MILEAGE

#include "TspBase.h"

enum EMInitType {
    EM_INIT_RANDOM,
    EM_INIT_EXTREMES,
    EM_INIT_FARTHEST_POINTS,
    EM_INIT_HULL // won't be working for now
};

Solution ExtraMileage(Instance *inst, enum EMInitType emType);

#endif // EXTRA_MILEAGE