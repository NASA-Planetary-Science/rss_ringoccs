#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_numerical.h>


/*  For floating type, set the error to 10^-8.                                */
#define EPS 1.0e-8F

float rssringoccs_Halleys_Method_Float(float x, float (*f)(float),
                                       float (*f_prime)(float),
                                       float (*f_2prime)(float),
                                       unsigned int max_iters)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    float dx, y, y_prime, y_2prime, denom;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    y = f(x);
    y_prime = f_prime(x);
    y_2prime = f_2prime(x);

    denom = 2.0F*y_prime*y_prime - y*y_2prime;

    /*  Check that the denominator is non-zero.                               */
    if (denom == 0.0F)
        return rssringoccs_NaN_F;

    /*  Compute the first iteration of Halley's method.                       */
    dx = 2.0F*y*y_prime/denom;
    x -= dx;

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_Float_Abs(dx) > EPS)
    {
        y = f(x);
        y_prime = f_prime(x);
        y_2prime = f_2prime(x);

        denom = 2*y_prime*y_prime - y*y_2prime;

        if (denom == 0.0F)
            return rssringoccs_NaN_F;

        dx = 2.0F*y*y_prime/denom;
        x -= dx;
        ++n;

        /*  Break if too many iterations have been run.                       */
        if (n > max_iters)
            break;
    }

    return x;
}

/*  Reset the error to 10^-16 for double and long double.                     */
#undef EPS
#define EPS 1.0e-16

double rssringoccs_Halleys_Method_Double(double x, double (*f)(double),
                                         double (*f_prime)(double),
                                         double (*f_2prime)(double),
                                         unsigned int max_iters)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double dx, y, y_prime, y_2prime, denom;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    y = f(x);
    y_prime = f_prime(x);
    y_2prime = f_2prime(x);

    denom = 2.0*y_prime*y_prime - y*y_2prime;

    /*  Check that the denominator is non-zero.                               */
    if (denom == 0.0)
        return rssringoccs_NaN;

    /*  Compute the first iteration of Halley's method.                       */
    dx = 2.0*y*y_prime/denom;
    x -= dx;

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_Double_Abs(dx) > EPS)
    {
        y = f(x);
        y_prime = f_prime(x);
        y_2prime = f_2prime(x);

        denom = 2.0*y_prime*y_prime - y*y_2prime;

        if (denom == 0.0)
            return rssringoccs_NaN;

        dx = 2.0*y*y_prime/denom;
        x -= dx;
        ++n;

        /*  Break if too many iterations have been run.                       */
        if (n > max_iters)
            break;
    }

    return x;
}

/*  Reset the error to 10^-16 for double and long double.                     */
#undef EPS
#define EPS 1.0e-16L

long double
rssringoccs_Halleys_Method_LDouble(long double x,
                                   long double (*f)(long double),
                                   long double (*f_prime)(long double),
                                   long double (*f_2prime)(long double),
                                   unsigned int max_iters)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    long double dx, y, y_prime, y_2prime, denom;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    y = f(x);
    y_prime = f_prime(x);
    y_2prime = f_2prime(x);

    denom = 2.0L*y_prime*y_prime - y*y_2prime;

    /*  Check that the denominator is non-zero.                               */
    if (denom == 0.0L)
        return rssringoccs_NaN_L;

    /*  Compute the first iteration of Halley's method.                       */
    dx = 2.0L*y*y_prime/denom;
    x -= dx;

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_LDouble_Abs(dx) > EPS)
    {
        y = f(x);
        y_prime = f_prime(x);
        y_2prime = f_2prime(x);

        denom = 2.0L*y_prime*y_prime - y*y_2prime;

        if (denom == 0.0L)
            return rssringoccs_NaN_L;

        dx = 2.0L*y*y_prime/denom;
        x -= dx;
        ++n;

        /*  Break if too many iterations have been run.                       */
        if (n > max_iters)
            break;
    }

    return x;
}

/*  Reset the error to 10^-16 for double and long double.                     */
#undef EPS
#define EPS 1.0e-16

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
