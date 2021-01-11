#include <stdlib.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  To avoid repeating the same code over and over again, define this macro   *
 *  to be used for all of the where_greater functions. Since the only thing   *
 *  that changes between the various functions is the type of the input       *
 *  pointer, the code is exactly the same.                                    */

#define CreateWhereFunction(type, Type, threshtype)                            \
unsigned long **rssringoccs_Where_Greater_##Type(type *data, unsigned long dim,\
                                        threshtype threshold)                  \
{                                                                              \
    /*  Declare necessary variables.                                         */\
    unsigned long i;                                                           \
    unsigned long *arr;                                                        \
    unsigned long **where;                                                     \
    unsigned long arr_size = 0;                                                \
                                                                               \
    /*  Create a dummy array to store the correct indices in for now.        */\
    arr = malloc(sizeof(*arr) * dim);                                          \
                                                                               \
    /*  Loop over the input array data to check which indices match.         */\
    for (i=0; i<dim; ++i)                                                      \
    {                                                                          \
        if (data[i] > threshold)                                               \
        {                                                                      \
            arr[arr_size] = i;                                                 \
            arr_size += 1;                                                     \
        }                                                                      \
    }                                                                          \
                                                                               \
    /*  Resize the array dummy to avoid hogging extra memory.                */\
    arr = realloc(arr, sizeof(*arr) * arr_size);                               \
                                                                               \
    /*  Declare a pointer to a pointer to containg all info to be returned.  */\
    where = malloc(sizeof(*where) * 2);                                        \
                                                                               \
    /*  The first index is the actual array of indices.                      */\
    where[0] = arr;                                                            \
                                                                               \
    /*  The second index is the size of the returning array.                 */\
    where[1] = malloc(sizeof(*where[1]));                                      \
    *where[1] = arr_size;                                                      \
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
CreateWhereFunction(long double, LDouble, long double)

/*  Erase the CreateWhereFunction macro.                                      */
#undef CreateWhereFunction

#define CreateWhereFunction(type, Type, threshtype)                            \
unsigned long **rssringoccs_Where_Lesser_##Type(type *data, unsigned long dim, \
                                       threshtype threshold)                   \
{                                                                              \
    /*  Declare necessary variables.                                         */\
    unsigned long i;                                                           \
    unsigned long *arr;                                                        \
    unsigned long **where;                                                     \
    unsigned long arr_size = 0;                                                \
                                                                               \
    /*  Create a dummy array to store the correct indices in for now.        */\
    arr = malloc(sizeof(*arr) * dim);                                          \
                                                                               \
    /*  Loop over the input array data to check which indices match.         */\
    for (i=0; i<dim; ++i)                                                      \
    {                                                                          \
        if (data[i] < threshold)                                               \
        {                                                                      \
            arr[arr_size] = i;                                                 \
            arr_size += 1;                                                     \
        }                                                                      \
    }                                                                          \
                                                                               \
    /*  Resize the array dummy to avoid hogging extra memory.                */\
    arr = realloc(arr, sizeof(*arr) * arr_size);                               \
                                                                               \
    /*  Declare a pointer to a pointer to containg all info to be returned.  */\
    where = malloc(sizeof(*where) * 2);                                        \
                                                                               \
    /*  The first index is the actual array of indices.                      */\
    where[0] = arr;                                                            \
                                                                               \
    /*  The second index is the size of the returning array.                 */\
    where[1] = malloc(sizeof(*where[1]));                                      \
    *where[1] = arr_size;                                                      \
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
CreateWhereFunction(long double, LDouble, long double)

#undef CreateWhereFunction

#define CreateWhereFunction(type, Type, threshtype)\
unsigned long **rssringoccs_Where_LesserGreater_##Type(type *data, unsigned long dim,\
                                              threshtype lower,                \
                                              threshtype upper)                \
{                                                                              \
    /*  Declare necessary variables.                                         */\
    unsigned long i;                                                           \
    unsigned long *arr;                                                        \
    unsigned long **where;                                                     \
    unsigned long arr_size = 0;                                                \
                                                                               \
    /*  Create a dummy array to store the correct indices in for now.        */\
    arr = malloc(sizeof(*arr) * dim);                                          \
                                                                               \
    /*  Loop over the input array data to check which indices match.         */\
    for (i=0; i<dim; ++i)                                                      \
    {                                                                          \
        if ((data[i] > lower) && (data[i] < upper))                            \
        {                                                                      \
            arr[arr_size] = i;                                                 \
            arr_size += 1;                                                     \
        }                                                                      \
    }                                                                          \
                                                                               \
    /*  Resize the array dummy to avoid hogging extra memory.                */\
    arr = realloc(arr, sizeof(*arr) * arr_size);                               \
                                                                               \
    /*  Declare a pointer to a pointer to containg all info to be returned.  */\
    where = malloc(sizeof(*where) * 2);                                        \
                                                                               \
    /*  The first index is the actual array of indices.                      */\
    where[0] = arr;                                                            \
                                                                               \
    /*  The second index is the size of the returning array.                 */\
    where[1] = malloc(sizeof(*where[1]));                                      \
    *where[1] = arr_size;                                                      \
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
CreateWhereFunction(long double, LDouble, long double)
