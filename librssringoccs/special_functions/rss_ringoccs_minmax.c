/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define MinFuncFromType(type, Type)                                            \
type rssringoccs_Min_##Type(type *arr, long n_elements)                        \
{                                                                              \
    type min = arr[0];                                                         \
    long i;                                                                    \
    for (i = 1; i < n_elements; ++i)                                           \
    {                                                                          \
        if (arr[i] < min)                                                      \
            min = arr[i];                                                      \
    }                                                                          \
                                                                               \
    return min;                                                                \
}

MinFuncFromType(float, Float)
MinFuncFromType(double, Double)
MinFuncFromType(long double, LDouble)
MinFuncFromType(char, Char)
MinFuncFromType(unsigned char, UChar)
MinFuncFromType(short, Short)
MinFuncFromType(unsigned short, UShort)
MinFuncFromType(int, Int)
MinFuncFromType(unsigned int, UInt)
MinFuncFromType(long, Long)
MinFuncFromType(unsigned long, ULong)

#define MaxFuncFromType(type, Type)                                            \
type rssringoccs_Max_##Type(type *arr, long n_elements)                        \
{                                                                              \
    type max = arr[0];                                                         \
    long i;                                                                    \
    for (i = 1; i < n_elements; ++i)                                           \
    {                                                                          \
        if (max < arr[i])                                                      \
            max = arr[i];                                                      \
    }                                                                          \
                                                                               \
    return max;                                                                \
}

MaxFuncFromType(float, Float)
MaxFuncFromType(double, Double)
MaxFuncFromType(long double, LDouble)
MaxFuncFromType(char, Char)
MaxFuncFromType(unsigned char, UChar)
MaxFuncFromType(short, Short)
MaxFuncFromType(unsigned short, UShort)
MaxFuncFromType(int, Int)
MaxFuncFromType(unsigned int, UInt)
MaxFuncFromType(long, Long)
MaxFuncFromType(unsigned long, ULong)
