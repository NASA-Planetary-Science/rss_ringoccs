/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_kaiser_bessel_modified(type, Type)                             \
type rssringoccs_Modified_Kaiser_Bessel_##Type(type x, type W, type alpha)     \
{                                                                              \
    type kbmd, arg, abs_x;                                                     \
                                                                               \
    abs_x = rssringoccs_Abs_##Type(x);                                         \
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
            arg = rssringoccs_Sqrt_##Type(1.0 - arg*arg);                      \
                                                                               \
            alpha *= ONE_PI;                                                   \
            kbmd = (rssringoccs_Bessel_I0_##Type(alpha * arg) - 1.0) /         \
                   (rssringoccs_Bessel_I0_##Type(alpha) - 1.0);                \
        }                                                                      \
        else                                                                   \
            kbmd = 0.0;                                                        \
    }                                                                          \
                                                                               \
    return kbmd;                                                               \
}

_define_kaiser_bessel_modified(float, Float)
_define_kaiser_bessel_modified(double, Double)
_define_kaiser_bessel_modified(long double, LongDouble)
