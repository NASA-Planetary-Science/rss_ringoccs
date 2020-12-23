/*  Needed for puts.                                                          */
#include <stdio.h>

/*  Needed for exit and NULL.                                                 */
#include <stdlib.h>

#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble
rssringoccs_CDouble_Poly_Deriv_Real_Coeffs(double *coeffs,
                                           unsigned int degree,
                                           unsigned int deriv,
                                           rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble poly;
    double temp, factor;
    unsigned int n;

    /*  If the input coefficient pointer is NULL, trying to access it will    *
     *  result in a segmentation fault. Check this and abort the computation  *
     *  if it's NULL.                                                         */
    if (coeffs == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tFunction: rssringoccs_ComplexPolyRealCoeffs\n\n"
             "The input coefficients pointer is NULL. Trying to access this\n"
             "will result in a segmentation fault. Aborting computation.\n\n");
        exit(0);
    }

    /*  Degree should be at least one, otherwise this is not a polynomial but *
     *  a constant. Check this. If degree is zero, we'll just return the      *
     *  zeroth coefficient (a constant polynomial).                           */
    if (degree < deriv)
        poly = rssringoccs_CDouble_Zero;
    else if (degree == deriv)
        poly = rssringoccs_CDouble_Rect(coeffs[degree], 0.0);
    else
    {
        /*  Set poly to a_{N}*z + a_{N-1} where N is the degree.              */
        factor = (double)rssringoccs_Falling_Factorial(degree, deriv);
        temp = factor*coeffs[degree];
        poly = rssringoccs_CDouble_Multiply_Real(temp, z);

        factor = (double)rssringoccs_Falling_Factorial(degree-1, deriv);
        temp   = factor*coeffs[degree-1];
        poly   = rssringoccs_CDouble_Add_Real(temp, poly);

        /*  Use Horner's method of polynomial computation.                    */
        for (n=2; n<=degree-deriv; ++n)
        {
            poly   =  rssringoccs_CDouble_Multiply(poly, z);
            factor = (double)rssringoccs_Falling_Factorial(degree-n, deriv);
            temp  *= coeffs[degree-n];
            poly   =  rssringoccs_CDouble_Add_Real(temp, poly);
        }
    }
    return poly;
}
