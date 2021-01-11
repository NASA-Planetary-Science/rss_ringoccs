/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_kaiser_bessel_modified(type, Type)                             \
type rssringoccs_##Type##_Modified_Kaiser_Bessel(type x, type W, type alpha)   \
{                                                                              \
    type kbmd, arg, abs_x;                                                     \
                                                                               \
    abs_x = rssringoccs_##Type##_Abs(x);                                       \
                                                                               \
    if (alpha == 0.0)                                                          \
    {                                                                          \
        if (abs_x < 0.5*W)                                                     \
            kbmd = 1.0;                                                        \
        else                                                                   \
            kbmd = 0.0;                                                        \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        if (abs_x < 0.5*W)                                                     \
        {                                                                      \
            arg = 2.0*abs_x/W;                                                 \
            arg = rssringoccs_##Type##_Sqrt(1.0 - arg*arg);                    \
                                                                               \
            alpha *= rssringoccs_One_Pi;                                                   \
            kbmd = (rssringoccs_##Type##_Bessel_I0(alpha * arg) - 1.0) /       \
                   (rssringoccs_##Type##_Bessel_I0(alpha) - 1.0);              \
        }                                                                      \
        else                                                                   \
            kbmd = 0.0;                                                        \
    }                                                                          \
                                                                               \
    return kbmd;                                                               \
}

_define_kaiser_bessel_modified(float, Float)
_define_kaiser_bessel_modified(double, Double)
_define_kaiser_bessel_modified(long double, LDouble)
