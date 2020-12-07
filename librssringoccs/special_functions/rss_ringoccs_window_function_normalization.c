/*  The C standard library header math.h is included here, as are aliases for *
 *  various functions, the macros INFINITY and NAN, as well as the max legal  *
 *  values for the exponential function which don't return INFINITY.          */
/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Complex variables and functions defined here.                             */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

#define _define_window_normalization(type, Type)                               \
type rssringoccs_##Type##_Window_Normalization(type *ker, long dim,            \
                                               type dx, type f_scale)          \
{                                                                              \
    /*  Declare variable for indexing.                                       */\
    long n;                                                                    \
    type out;                                                                  \
                                                                               \
    /*  Compute the Free-Space integral.                                     */\
    double T1 = 0.0;                                                           \
                                                                               \
    for (n=0; n<dim; ++n)                                                      \
        T1 += ker[n];                                                          \
                                                                               \
    T1 = rssringoccs_##Type##_Abs(T1 * dx);                                    \
                                                                               \
    /* Return the normalization factor.                                      */\
    out = rssringoccs_Sqrt_Two * f_scale / T1;                                               \
    return out;                                                                \
}

_define_window_normalization(float, Float)
_define_window_normalization(double, Double)
_define_window_normalization(long double, LDouble)


double
rssringoccs_Complex_Window_Normalization(rssringoccs_ComplexDouble *ker,
                                         long dim, double dx, double f_scale)
{
    /*  Declare variable for indexing.                                        */
    long n;
    double abs_T1, out;

    /*  Compute the Free-Space integral.                                      */
    rssringoccs_ComplexDouble T1 = rssringoccs_CDouble_Zero;

    for (n=0; n<dim; ++n)
        T1 = rssringoccs_CDouble_Add(T1, ker[n]);

    abs_T1 = dx*rssringoccs_CDouble_Abs(T1);

    /* Return the normalization factor.                                       */
    out = rssringoccs_Sqrt_Two * f_scale / abs_T1;
    return out;
}
