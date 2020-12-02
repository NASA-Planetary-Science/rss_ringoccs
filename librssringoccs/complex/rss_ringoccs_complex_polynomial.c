/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************
 *                     rss_ringoccs_complex_polynomial                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Efficiently compute polynomials of a complex variable using Horner's  *
 *      method. Two functions are provided allowing for real coefficients and *
 *      complex coefficients.                                                 *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 *  1.) stdio.h:                                                              *
 *          The standard C library for input and output. This is needed for   *
 *          the "puts" function for printing error messages should the        *
 *          program fail.                                                     *
 *  2.) stdlib:                                                               *
 *          The C standard library header file. This contains the "exit"      *
 *          function which is used to terminate a program if an error occurs. *
 *  2.) rss_ringoccs_complex.h:                                               *
 *          Header file where rssringoccs_ComplexDouble is defined, as well   *
 *          as the prototypes for the functions.                              *
 ******************************************************************************
 *                                 WARNINGS                                   *
 *  1.) These functions use "exit" to abort a computation if an error occurs. *
 *      This is to prevent the possibility of a segmentation fault, which     *
 *      often prints out a less-than-helpful error message. Improper use of   *
 *      this function will cause your program to crash (However, the error    *
 *      message will tell you why).                                           *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 12, 2020                                             *
 ******************************************************************************/

/*  Needed for the 'puts' function.                                           */
#include <stdio.h>

/*  And the 'exit' function is located here.                                  */
#include <stdlib.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Use Horner's method to compute a polynomial of a complex variable z with  *
 *  real coefficients.                                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Poly_Real_Coeffs(double *coeffs, unsigned int degree,
                                           rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble poly;
    unsigned int n;

    /*  If the input coefficient pointer is NULL, trying to access it will    *
     *  result in a segmentation fault. Check this and abort the computation  *
     *  if it's NULL.                                                         */
    if (coeffs == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tFunction: rssringoccs_Complex_Poly_Real_Coeffs\n\n"
             "The input coefficients pointer is NULL. Trying to access this\n"
             "will result in a segmentation fault. Aborting computation.\n\n");
        exit(0);
    }

    /*  Degree should be at least one, otherwise this is not a polynomial but *
     *  a constant. Check this. If degree is zero, we'll just return the      *
     *  zeroth coefficient (a constant polynomial).                           */
    if (degree == 0)
    {
        poly = rssringoccs_CDouble_Rect(coeffs[0], 0.0);
        return poly;
    }

    /*  Set poly equal to a_{N}*z, where N is the degree of the polynomial.   */
    poly = rssringoccs_CDouble_Multiply_Real(coeffs[degree], z);

    /*  Reset poly to a_{N}*z + a_{N-1}.                                      */
    poly = rssringoccs_CDouble_Add_Real(coeffs[degree-1], poly);

    /*  Use Horner's method of polynomial computation.                        */
    for (n=2; n<=degree; ++n)
    {
        /*  Use Horner's method with the current complex coefficients.        */
        poly = rssringoccs_CDouble_Multiply(z, poly);
        poly = rssringoccs_CDouble_Add_Real(coeffs[degree-n], poly);
    }

    return poly;
}

/*  Use Horner's method to compute a polynomial of a complex variable z with  *
 *  complex coefficients.                                                     */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Poly_Complex_Coeffs(rssringoccs_ComplexDouble *coeffs,
                                              unsigned int degree,
                                              rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble poly;
    unsigned int n;

    /*  If the input coefficient pointer is NULL, trying to access it will    *
     *  result in a segmentation fault. Check this and abort the computation  *
     *  if it's NULL.                                                         */
    if (coeffs == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tFunction: rssringoccs_Complex_Poly_Complex_Coeffs\n\n"
             "The input coefficients pointer is NULL. Trying to access this\n"
             "will result in a segmentation fault. Aborting computation.\n\n");
        exit(0);
    }

    /*  Degree should be at least one, otherwise this is not a polynomial but *
     *  a constant. Check this. If degree is zero, we'll just return the      *
     *  zeroth coefficient (a constant polynomial).                           */
    if (degree == 0)
        return coeffs[0];

    /*  Set poly equal to a_{N}*z, where N is the degree of the polynomial.   */
    poly = rssringoccs_CDouble_Multiply(coeffs[degree], z);

    /*  Reset poly to a_{N}*z + a_{N-1}.                                      */
    poly = rssringoccs_CDouble_Add(poly, coeffs[degree-1]);

    /*  Use Horner's method of polynomial computation.                        */
    for (n=2; n<=degree; ++n)
    {
        poly = rssringoccs_CDouble_Multiply(z, poly);
        poly = rssringoccs_CDouble_Add(poly, coeffs[degree-n]);
    }

    return poly;
}
