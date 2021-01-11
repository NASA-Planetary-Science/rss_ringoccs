/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_resolution_inverse(type, Type)                                 \
type rssringoccs_##Type##_Resolution_Inverse(type x)                           \
{                                                                              \
    type P1, P2, out;                                                          \
                                                                               \
    if (x <= 1.0)                                                              \
        out = rssringoccs_NaN;                                                 \
    else if (x < rssringoccs_Infinity)                                         \
    {                                                                          \
        P1 = x/(1.0-x);                                                        \
        P2 = P1*rssringoccs_##Type##_Exp(P1);                                  \
        out = rssringoccs_##Type##_LambertW(P2)-P1;                            \
    }                                                                          \
    else                                                                       \
        out = 0.0;                                                             \
                                                                               \
    return out;                                                                \
}

_define_resolution_inverse(float, Float)
_define_resolution_inverse(double, Double)
_define_resolution_inverse(long double, LDouble)
