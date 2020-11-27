/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_sinc(type, Type)                                               \
type rssringoccs_##Type##_Sinc(type x)                                         \
{                                                                              \
    type y;                                                                    \
                                                                               \
    if (x == 0.0)                                                              \
        y = 1.0;                                                               \
    else                                                                       \
        y = rssringoccs_##Type##_Sin(x)/x;                                     \
                                                                               \
    return y;                                                                  \
}

_define_sinc(float, Float)
_define_sinc(double, Double)
_define_sinc(long double, LongDouble)
