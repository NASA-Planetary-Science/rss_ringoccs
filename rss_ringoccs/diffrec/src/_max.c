#include "special_functions.h"

#ifdef MaxFuncFromType
#undef MaxFuncFromType
#endif

#define MaxFuncFromType(type, Type)\
type Max_##Type(type *arr, long n_elements)\
{\
    type max = arr[0];\
    long i;\
    for (i = 1; i < n_elements; ++i)\
    {\
        if (arr[i] > max) max = arr[i];\
    }\
    return max;\
}

MaxFuncFromType(float, Float);
MaxFuncFromType(double, Double);
MaxFuncFromType(long double, Long_Double);
MaxFuncFromType(short, Short);
MaxFuncFromType(int, Int);
MaxFuncFromType(long, Long);
MaxFuncFromType(long long, Long_Long);
