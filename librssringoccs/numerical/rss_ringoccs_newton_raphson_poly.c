#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_numerical.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>

#define EPS 1.0e-16

rssringoccs_ComplexDouble
rssringoccs_Newton_Raphson_CDouble_Poly_Real(
    rssringoccs_ComplexDouble z, double *coeffs, unsigned int degree,
    unsigned int max_iters
)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    rssringoccs_ComplexDouble dz, w, wp;
    rssringoccs_Bool comp;
    unsigned int n;

    /*  Evaluate the perturbation term, then compute the next iteration.      */
    w = rssringoccs_CDouble_Poly_Real_Coeffs(coeffs, degree, z);
    wp = rssringoccs_CDouble_Poly_Deriv_Real_Coeffs(coeffs, degree, 1, z);

    /*  If the derivative is zero at your initial guess, Newton-Raphson       *
     *  fails. Return Not-a-Number in this case.                              */
    comp = rssringoccs_CDouble_Compare(wp, rssringoccs_CDouble_Zero);
    if (comp)
        return rssringoccs_CDouble_NaN;

    /*  Compute the first iteration of Newton-Raphson.                        */
    dz = rssringoccs_CDouble_Divide(w, wp);
    z  = rssringoccs_CDouble_Subtract(z, dz);

    /*  The first iteration has been computed above, so set n to 1.           */
    n = 1;

    /*  Continuing this computation until the error is below the threshold.   */
    while(rssringoccs_CDouble_Abs(dz) > EPS)
    {
        w = rssringoccs_CDouble_Poly_Real_Coeffs(coeffs, degree, z);
        wp = rssringoccs_CDouble_Poly_Deriv_Real_Coeffs(coeffs, degree,
                                                              1, z);

        comp = rssringoccs_CDouble_Compare(wp, rssringoccs_CDouble_Zero);

        if (comp)
            return rssringoccs_CDouble_NaN;

        dz = rssringoccs_CDouble_Divide(w, wp);
        z  = rssringoccs_CDouble_Subtract(z, dz);
        ++n;

        /*  Break if too many iterations have been run.                       */
        if (n > max_iters)
            break;
    }

    return z;
}

#undef EPS
