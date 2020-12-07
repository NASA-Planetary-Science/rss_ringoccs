#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_numerical.h>

/*  For floating type, set the error to 10^-8.                                */
#define EPS 1.0e-8F

float rssringoccs_Newton_Raphson_Float(float x, float (*f)(float),
                                       float (*f_prime)(float),
                                       unsigned int max_iters)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    float dx, y, y_prime;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    y = (*f)(x);
    y_prime = (*f_prime)(x);

    /*  If the derivative is zero at your initial guess, Newton-Raphson       *
     *  fails. Return Not-a-Number in this case.                              */
    if (y_prime == 0.0F)
        return rssringoccs_NaN_F;

    /*  Compute the first iteration of Newton-Raphson.                        */
    dx = y/y_prime;
    x -= dx;

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_Float_Abs(dx) > EPS)
    {
        y = (*f)(x);
        y_prime = (*f_prime)(x);

        if (y_prime == 0.0F)
            return rssringoccs_NaN_F;

        dx = y/y_prime;
        x -= dx;
        ++n;

        /*  Break if too many iterations have been run.                      */
        if (n > max_iters)
            break;
    }

    return x;
}

/*  Reset the error to 10^-16 for double and long double.                     */
#undef EPS
#define EPS 1.0e-16

double rssringoccs_Newton_Raphson_Double(double x, double (*f)(double),
                                         double (*f_prime)(double),
                                         unsigned int max_iters)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double dx, y, y_prime;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    y = (*f)(x);
    y_prime = (*f_prime)(x);

    /*  If the derivative is zero at your initial guess, Newton-Raphson       *
     *  fails. Return Not-a-Number in this case.                              */
    if (y_prime == 0.0)
        return rssringoccs_NaN;

    /*  Compute the first iteration of Newton-Raphson.                        */
    dx = y/y_prime;
    x -= dx;

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_Double_Abs(dx) > EPS)
    {
        y = (*f)(x);
        y_prime = (*f_prime)(x);

        if (y_prime == 0.0)
            return rssringoccs_NaN;

        dx = y/y_prime;
        x -= dx;
        ++n;

        /*  Break if too many iterations have been run.                      */
        if (n > max_iters)
            break;
    }

    return x;
}

rssringoccs_ComplexDouble
rssringoccs_Newton_Raphson_Complex(
    rssringoccs_ComplexDouble z,
    rssringoccs_ComplexDouble (*f)(rssringoccs_ComplexDouble),
    rssringoccs_ComplexDouble (*f_prime)(rssringoccs_ComplexDouble),
    unsigned int max_iters
)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    rssringoccs_ComplexDouble dz, w, w_prime;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    w = (*f)(z);
    w_prime = (*f_prime)(z);

    /*  If the derivative is zero at your initial guess, Newton-Raphson       *
     *  fails. Return Not-a-Number in this case.                              */
    if (rssringoccs_CDouble_Compare(w_prime, rssringoccs_CDouble_Zero))
        return rssringoccs_CDouble_NaN;

    /*  Compute the first iteration of Newton-Raphson.                        */
    dz = rssringoccs_CDouble_Divide(w, w_prime);
    z  = rssringoccs_CDouble_Subtract(z, dz);

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_CDouble_Abs(dz) > EPS)
    {
        w = (*f)(z);
        w_prime = (*f_prime)(z);

        if (rssringoccs_CDouble_Compare(w_prime, rssringoccs_CDouble_Zero))
            return rssringoccs_CDouble_NaN;

        dz = rssringoccs_CDouble_Divide(w, w_prime);
        z  = rssringoccs_CDouble_Subtract(z, dz);
        ++n;

        /*  Break if too many iterations have been run.                       */
        if (n > max_iters)
            break;
    }

    return z;
}
