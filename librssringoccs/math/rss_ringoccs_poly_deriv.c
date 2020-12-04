/*  Needed for puts.                                                          */
#include <stdio.h>

/*  Needed for exit and NULL.                                                 */
#include <stdlib.h>

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

float rssringoccs_Real_Poly_Deriv_Float_Coeffs(float *coeffs,
                                               unsigned int degree,
                                               unsigned int deriv, float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float poly;
    unsigned int n, N;

    N = degree - deriv;

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

    /*  Set poly to a_{N}*z + a_{N-1} where N is the degree.                  */
    poly = rssringoccs_Falling_Factorial(degree, deriv)*coeffs[0];

    /*  Use Horner's method of polynomial computation.                        */
    for (n=1; n<=N; ++n)
        poly = x*poly +
            rssringoccs_Falling_Factorial(degree-n, deriv)*coeffs[degree-n];

    return poly;
}

double rssringoccs_Real_Poly_Deriv_Double_Coeffs(double *coeffs,
                                                 unsigned int degree,
                                                 unsigned int deriv, double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double poly;
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
    if (degree == 0)
        poly = coeffs[0];
    else

    /*  Set poly to a_{N}*z + a_{N-1} where N is the degree.                  */
    poly = rssringoccs_Falling_Factorial(degree, deriv)*coeffs[0];

    /*  Use Horner's method of polynomial computation.                        */
    for (n=1; n<=degree; ++n)
        poly = x*poly +
            rssringoccs_Falling_Factorial(degree-n, deriv)*coeffs[degree-n];

    return poly;
}

long double
rssringoccs_Real_Poly_Deriv_LDouble_Coeffs(long double *coeffs,
                                           unsigned int degree,
                                           unsigned int deriv,
                                           long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double poly;
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

    /*  Set poly to a_{N}*z + a_{N-1} where N is the degree.                  */
    poly = rssringoccs_Falling_Factorial(degree, deriv)*coeffs[0];

    /*  Use Horner's method of polynomial computation.                        */
    for (n=1; n<=degree; ++n)
        poly = x*poly +
            rssringoccs_Falling_Factorial(degree-n, deriv)*coeffs[degree-n];

    return poly;
}
