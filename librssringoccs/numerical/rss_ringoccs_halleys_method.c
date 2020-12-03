#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_numerical.h>

#define __do_halleys_method(type, Type)                                        \
type rssringoccs_Halleys_Method_##Type(type x, type (*f)(type),                \
                                       type (*f_prime)(type),                  \
                                       type (*f_2prime)(type),                 \
                                       unsigned int max_iters)                 \
{                                                                              \
    /*  Declare necessary variables. C89 requires this at the top.           */\
    type dx, y, y_prime, y_2prime, denom;                                      \
    unsigned int n;                                                            \
                                                                               \
    /*  Evaluate the perturbation term, then compute the next iteration.     */\
    y = f(x);                                                                  \
    y_prime = f_prime(x);                                                      \
    y_2prime = f_2prime(x);                                                    \
                                                                               \
    denom = 2*y_prime*y_prime - y*y_2prime;                                    \
                                                                               \
    /*  Check that the denominator is non-zero.                              */\
    if (denom == 0.0)                                                          \
        return rssringoccs_NaN;                                                \
                                                                               \
    /*  Compute the first iteration of Halley's method.                      */\
    dx = 2*y*y_prime/denom;                                                    \
    x -= dx;                                                                   \
                                                                               \
    /*  The first iteration has been computed above, so set n to 1.          */\
    n = 1;                                                                     \
                                                                               \
    /*  Continuing this computation until the error is below the threshold.  */\
    while(rssringoccs_##Type##_Abs(dx) > EPS)                                  \
    {                                                                          \
        y = f(x);                                                              \
        y_prime = f_prime(x);                                                  \
        y_2prime = f_2prime(x);                                                \
                                                                               \
        denom = 2*y_prime*y_prime - y*y_2prime;                                \
                                                                               \
        if (denom == 0.0)                                                      \
            return rssringoccs_NaN;                                            \
                                                                               \
        dx = 2*y*y_prime/denom;                                                \
        x -= dx;                                                               \
        ++n;                                                                   \
                                                                               \
        /*  Break if too many iterations have been run.                      */\
        if (n > max_iters)                                                     \
            break;                                                             \
    }                                                                          \
                                                                               \
    return x;                                                                  \
}

/*  For floating type, set the error to 10^-8.                                */
#define EPS 1.0e-8

__do_halleys_method(float, Float)

/*  Reset the error to 10^-16 for double and long double.                     */
#undef EPS
#define EPS 1.0e-16

__do_halleys_method(double, Double)
__do_halleys_method(long double, LDouble)

rssringoccs_ComplexDouble
rssringoccs_Halleys_Method_Complex(
    rssringoccs_ComplexDouble z,
    rssringoccs_ComplexDouble (*f)(rssringoccs_ComplexDouble),
    rssringoccs_ComplexDouble (*f_prime)(rssringoccs_ComplexDouble),
    rssringoccs_ComplexDouble (*f_2prime)(rssringoccs_ComplexDouble),
    unsigned int max_iters
)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    rssringoccs_ComplexDouble dz, w, w_prime, w_2prime, denom;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    w = f(z);
    w_prime = f_prime(z);
    w_2prime = f_2prime(z);

    denom = rssringoccs_CDouble_Subtract(
        rssringoccs_CDouble_Multiply_Real(
            2.0 , rssringoccs_CDouble_Multiply(w_prime, w_prime)
        ),
        rssringoccs_CDouble_Multiply(w, w_2prime)
    );

    /*  Check that the denominator is non-zero.                               */
    if (rssringoccs_CDouble_Compare(denom, rssringoccs_CDouble_Zero))
        return rssringoccs_CDouble_NaN;

    /*  Compute the first iteration of Newton-Raphson.                        */
    dz = rssringoccs_CDouble_Divide(
        rssringoccs_CDouble_Multiply_Real(
            2.0, rssringoccs_CDouble_Multiply(w, w_prime)
        ),
        denom
    );

    z = rssringoccs_CDouble_Subtract(z, dz);

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_CDouble_Abs(dz) > EPS)
    {
        w = f(z);
        w_prime = f_prime(z);
        w_2prime = f_2prime(z);

        denom = rssringoccs_CDouble_Subtract(
            rssringoccs_CDouble_Multiply_Real(
                2.0 , rssringoccs_CDouble_Multiply(w_prime, w_prime)
            ),
            rssringoccs_CDouble_Multiply(w, w_2prime)
        );

        /*  Check that the denominator is non-zero.                           */
        if (rssringoccs_CDouble_Compare(denom, rssringoccs_CDouble_Zero))
            return rssringoccs_CDouble_NaN;

        dz = rssringoccs_CDouble_Divide(
            rssringoccs_CDouble_Multiply_Real(
                2.0, rssringoccs_CDouble_Multiply(w, w_prime)
            ),
            denom
        );

        z = rssringoccs_CDouble_Subtract(z, dz);
        ++n;

        /*  Break if too many iterations have been run.                       */
        if (n > max_iters)
            break;
    }

    return z;
}
