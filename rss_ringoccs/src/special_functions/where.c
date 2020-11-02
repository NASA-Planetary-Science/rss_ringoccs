#include <stdlib.h>

/*  Prototypes for these functions declared here.                             */
#include "rss_ringoccs_special_functions.h"

/*  To avoid repeating the same code over and over again, define this macro   *
 *  to be used for all of the where_greater functions. Since the only thing   *
 *  that changes between the various functions is the type of the input       *
 *  pointer, the code is exactly the same.                                    */

#define CreateWhereFunction(type, Type, threshtype)                            \
long **rssringoccs_Where_Greater_##Type(type *data, long dim,                  \
                                        threshtype threshold)                  \
{                                                                              \
    /*  Declare necessary variables.                                         */\
    long i;                                                                    \
    long *dummy;                                                               \
    long *arr;                                                                 \
    long **where;                                                              \
    long arr_size = 0;                                                         \
                                                                               \
    /*  Create a dummy array to store the correct indices in for now.        */\
    dummy = (long *)malloc(sizeof(long) * dim);                                \
                                                                               \
    /*  Loop over the input array data to check which indices match.         */\
    for (i=0; i<dim; ++i)                                                      \
    {                                                                          \
        if (data[i] > threshold)                                               \
        {                                                                      \
            dummy[arr_size] = i;                                               \
            arr_size += 1;                                                     \
        }                                                                      \
    }                                                                          \
                                                                               \
    /*  Resize the array dummy to avoid hogging extra memory.                */\
    arr = (long *)malloc(sizeof(long) * arr_size);                             \
    for (i=0; i<arr_size; ++i)                                                 \
        arr[i] = dummy[i];                                                     \
                                                                               \
    /*  Since arr now contains all of the necessary info, free dummy.        */\
    free(dummy);                                                               \
                                                                               \
    /*  Declare a pointer to a pointer to containg all info to be returned.  */\
    where = malloc(sizeof(*where) * 2);                                        \
                                                                               \
    /*  The first index is the actual array of indices.                      */\
    where[0] = arr;                                                            \
                                                                               \
    /*  The second index is the size of the returning array.                 */\
    where[1] = &arr_size;                                                      \
    return where;                                                              \
}
 
CreateWhereFunction(char, Char, double)
CreateWhereFunction(unsigned char, UChar, double)
CreateWhereFunction(short, Short, double)
CreateWhereFunction(unsigned short, UShort, double)
CreateWhereFunction(int, Int, double)
CreateWhereFunction(unsigned int, UInt, double)
CreateWhereFunction(long, Long, double)
CreateWhereFunction(unsigned long, ULong, double)
CreateWhereFunction(float, Float, float)
CreateWhereFunction(double, Double, double)
CreateWhereFunction(long double, LongDouble, long double)

/*  Erase the CreateWhereFunction macro.                                      */
#undef CreateWhereFunction

#define CreateWhereFunction(type, Type, threshtype)                            \
long **rssringoccs_Where_Lesser_##Type(type *data, long dim,                   \
                                       threshtype threshold)                   \
{                                                                              \
    /*  Declare necessary variables.                                         */\
    long i;                                                                    \
    long *dummy;                                                               \
    long *arr;                                                                 \
    long **where;                                                              \
    long arr_size = 0;                                                         \
                                                                               \
    /*  Create a dummy array to store the correct indices in for now.        */\
    dummy = (long *)malloc(sizeof(long) * dim);                                \
                                                                               \
    /*  Loop over the input array data to check which indices match.         */\
    for (i=0; i<dim; ++i)                                                      \
    {                                                                          \
        if (data[i] < threshold){                                              \
            dummy[arr_size] = i;                                               \
            arr_size += 1;                                                     \
        }                                                                      \
    }                                                                          \
                                                                               \
    /*  Resize the array dummy to avoid hogging extra memory.                */\
    arr = (long *)malloc(sizeof(long) * arr_size);                             \
    for (i=0; i<arr_size; ++i)                                                 \
        arr[i] = dummy[i];                                                     \
                                                                               \
    /*  Since arr now contains all of the necessary info, free dummy.        */\
    free(dummy);                                                               \
                                                                               \
    /*  Declare a pointer to a pointer to containg all info to be returned.  */\
    where = (long **)malloc(sizeof(long) * 2);                                 \
                                                                               \
    /*  The first index is the actual array of indices.                      */\
    where[0] = arr;                                                            \
                                                                               \
    /*  The second index is the size of the returning array.                 */\
    where[1] = &arr_size;                                                      \
    return where;                                                              \
}

CreateWhereFunction(char, Char, double)
CreateWhereFunction(unsigned char, UChar, double)
CreateWhereFunction(short, Short, double)
CreateWhereFunction(unsigned short, UShort, double)
CreateWhereFunction(int, Int, double)
CreateWhereFunction(unsigned int, UInt, double)
CreateWhereFunction(long, Long, double)
CreateWhereFunction(unsigned long, ULong, double)
CreateWhereFunction(float, Float, float)
CreateWhereFunction(double, Double, double)
CreateWhereFunction(long double, LongDouble, long double)

#undef CreateWhereFunction

#define CreateWhereFunction(type, Type, threshtype)\
long **rssringoccs_Where_LesserGreater_##Type(type *data, long dim,            \
                                              threshtype lower,                \
                                              threshtype upper)                \
{                                                                              \
    /*  Declare necessary variables.                                         */\
    long i;                                                                    \
    long *dummy;                                                               \
    long *arr;                                                                 \
    long **where;                                                              \
    long arr_size = 0;                                                         \
                                                                               \
    /*  Create a dummy array to store the correct indices in for now.        */\
    dummy = (long *)malloc(sizeof(long) * dim);                                \
                                                                               \
    /*  Loop over the input array data to check which indices match.         */\
    for (i=0; i<dim; ++i) {                                                    \
        if ((data[i] > lower) && (data[i] < upper))                            \
        {                                                                      \
            dummy[arr_size] = i;                                               \
            arr_size += 1;                                                     \
        }                                                                      \
    }                                                                          \
                                                                               \
    /*  Resize the array dummy to avoid hogging extra memory.                */\
    arr = (long *)malloc(sizeof(long) * arr_size);                             \
    for (i=0; i<arr_size; ++i)                                                 \
        arr[i] = dummy[i];                                                     \
                                                                               \
    /*  Since arr now contains all of the necessary info, free dummy.        */\
    free(dummy);                                                               \
                                                                               \
    /*  Declare a pointer to a pointer to containg all info to be returned.  */\
    where = (long **)malloc(sizeof(long) * 2);                                 \
                                                                               \
    /*  The first index is the actual array of indices.                      */\
    where[0] = arr;                                                            \
                                                                               \
    /*  The second index is the size of the returning array.                 */\
    where[1] = &arr_size;                                                      \
    return where;                                                              \
}

CreateWhereFunction(char, Char, double)
CreateWhereFunction(unsigned char, UChar, double)
CreateWhereFunction(short, Short, double)
CreateWhereFunction(unsigned short, UShort, double)
CreateWhereFunction(int, Int, double)
CreateWhereFunction(unsigned int, UInt, double)
CreateWhereFunction(long, Long, double)
CreateWhereFunction(unsigned long, ULong, double)
CreateWhereFunction(float, Float, float)
CreateWhereFunction(double, Double, double)
CreateWhereFunction(long double, LongDouble, long double)
