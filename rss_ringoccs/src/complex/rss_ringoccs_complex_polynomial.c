/*  Needed for the 'puts' function.                                           */
#include <stdio.h>

/*  And the 'exit' function is located here.                                  */
#include <stdlib.h>


/*  Where the prototypes are declared and where complex types are defined.    */
#include "rss_ringoccs_complex.h"

/*  Use Horner's method to compute a polynomial of a complex variable z with  *
 *  real coefficients.                                                        */
rssringoccs_ComplexDouble
rssringoccs_Complex_Poly_Real_Coeffs(double *coeffs, unsigned int degree,
                                     rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble poly, next_coeff;
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
    {
        poly = rssringoccs_Complex_Rect(coeffs[0], 0.0);
        return poly;
    }

    /*  Compute the first term in Horner's method. We need the coefficient to *
     *  be a complex number so we can add it to the complex z value.          */
    next_coeff = rssringoccs_Complex_Rect(coeffs[degree-1], 0.0);

    /*  Set poly equal to a_{N}*z, whre N is the degree of the polynomial.    */
    poly = rssringoccs_Complex_Scale(coeffs[degree], z);

    /*  Reset poly to a_{N}*z + a_{N-1}.                                      */
    poly = rssringoccs_Complex_Add(poly, next_coeff);

    /*  Use Horner's method of polynomial computation.                        */
    for (n=2; n<=degree; ++n)
    {
        /*  Convert the current coefficient to a complex value.               */
        next_coeff = rssringoccs_Complex_Rect(coeffs[degree-n], 0.0);

        /*  Use Horner's method with the current complex coefficients.        */
        poly = rssringoccs_Complex_Multiply(z, poly);
        poly = rssringoccs_Complex_Add(poly, next_coeff);
    }

    return poly;
}
