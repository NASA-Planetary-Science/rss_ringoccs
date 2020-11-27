/*  The C standard library header math.h is included here, as are aliases for *
 *  various functions, the macros INFINITY and NAN, as well as the max legal  *
 *  values for the exponential function which don't return INFINITY.          */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Complex variables and functions defined here.                             */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_coss_window(type, Type)                                        \
type rssringoccs_Coss_Window_##Type(type x, type W)                            \
{                                                                              \
    /*  Declare necessary variables. C89 requires declaring these at the top.*/\
    type abs_x, coss_x;                                                        \
                                                                               \
    /*  Get the absolute value of x.                                         */\
    abs_x = rssringoccs_Abs_##Type(x);                                         \
                                                                               \
    /*  Compute the cosine squared window function.                          */\
    if (abs_x <= W/2.0)                                                        \
    {                                                                          \
        abs_x *= ONE_PI/W;                                                     \
        coss_x = rssringoccs_##Type##_Cos(abs_x);                              \
        coss_x *= coss_x;                                                      \
    }                                                                          \
    else                                                                       \
        coss_x = 0.0;                                                          \
                                                                               \
    return coss_x;                                                             \
}

_define_coss_window(float, Float)
_define_coss_window(double, Double)
_define_coss_window(long double, LongDouble)
