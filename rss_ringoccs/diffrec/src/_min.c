#include "special_functions.h"

#ifdef MinFuncFromType
#undef MinFuncFromType
#endif

#define MinFuncFromType(type, Type)\
type Min_##Type(type *arr, long n_elements)\
{\
    type max = arr[0];\
    long i;\
    for (i = 1; i < n_elements; ++i)\
    {\
        if (arr[i] < max) max = arr[i];\
    }\
    return max;\
}

MinFuncFromType(float, Float);
MinFuncFromType(double, Double);
MinFuncFromType(long double, Long_Double);
MinFuncFromType(short, Short);
MinFuncFromType(int, Int);
MinFuncFromType(long, Long);
MinFuncFromType(long long, Long_Long);
