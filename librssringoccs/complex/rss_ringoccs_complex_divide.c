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
 *                      rss_ringoccs_complex_divide                           *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex division.                        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 30, 2020                                             *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  In C99, since _Complex is a built-in data type, given double _Complex z1  *
 *  and double _Complex z2, you can just do z1 / z2. Structs cannot be        *
 *  divided so we need a function for computing this.                         */

/*  Single precision complex division.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Divide(rssringoccs_ComplexFloat z0,
                                rssringoccs_ComplexFloat z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat div;
    float real0, imag0, real1, imag1, denom;
    float div_re, div_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_CFloat_Real_Part(z0);
    real1 = rssringoccs_CFloat_Real_Part(z1);
    imag0 = rssringoccs_CFloat_Imag_Part(z0);
    imag1 = rssringoccs_CFloat_Imag_Part(z1);

    /*  The denominator for both real and imaginary parts is |z1|^2.          */
    denom = 1.0 / rssringoccs_CFloat_Abs_Squared(z1);

    /*  We compute based on the fact that z0/z1 = z0 * (z1)^-1 and use the    *
     *  formular for the reciprocal of a complex number.                      */
    div_re = (real0*real1 + imag0*imag1)*denom;
    div_im = (imag0*real1 - real0*imag1)*denom;

    /*  Use rssringoccs_CFloat_Rect to create the output and return.    */
    div = rssringoccs_CFloat_Rect(div_re, div_im);
    return div;
}
/*  End of rssringoccs_CFloat_Divide.                                   */

/*  Double precision complex division.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Divide(rssringoccs_ComplexDouble z0,
                                 rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble div;
    double real0, imag0, real1, imag1, denom;
    double div_re, div_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_CDouble_Real_Part(z0);
    real1 = rssringoccs_CDouble_Real_Part(z1);
    imag0 = rssringoccs_CDouble_Imag_Part(z0);
    imag1 = rssringoccs_CDouble_Imag_Part(z1);

    /*  The denominator for both real and imaginary parts is |z1|^2.          */
    denom = 1.0 / rssringoccs_CDouble_Abs_Squared(z1);

    /*  We compute based on the fact that z0/z1 = z0 * (z1)^-1 and use the    *
     *  formular for the reciprocal of a complex number.                      */
    div_re = (real0*real1 + imag0*imag1)*denom;
    div_im = (imag0*real1 - real0*imag1)*denom;

    /*  Use rssringoccs_CDouble_Rect to create the output and return.   */
    div = rssringoccs_CDouble_Rect(div_re, div_im);
    return div;
}
/*  End of rssringoccs_CDouble_Divide.                                  */

/*  Long double precision complex division.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Divide(rssringoccs_ComplexLongDouble z0,
                                     rssringoccs_ComplexLongDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble div;
    long double real0, imag0, real1, imag1, denom;
    long double div_re, div_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_CLDouble_Real_Part(z0);
    real1 = rssringoccs_CLDouble_Real_Part(z1);
    imag0 = rssringoccs_CLDouble_Imag_Part(z0);
    imag1 = rssringoccs_CLDouble_Imag_Part(z1);

    /*  The denominator for both real and imaginary parts is |z1|^2.          */
    denom = 1.0 / rssringoccs_CLDouble_Abs_Squared(z1);

    /*  We compute based on the fact that z0/z1 = z0 * (z1)^-1 and use the    *
     *  formular for the reciprocal of a complex number.                      */
    div_re = (real0*real1 + imag0*imag1)*denom;
    div_im = (imag0*real1 - real0*imag1)*denom;

    /*  Use rssringoccs_CLDouble_Rect and return output.             */
    div = rssringoccs_CLDouble_Rect(div_re, div_im);
    return div;
}
/*  End of rssringoccs_CLDouble_Divide.                              */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here we have complex.h support so we'll use the / symbol.       */

/*  Single precision complex division.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Divide(rssringoccs_ComplexFloat z0,
                                rssringoccs_ComplexFloat z1)
{
    return z0/z1;
}
/*  End of rssringoccs_CFloat_Divide.                                   */

/*  Double precision complex division.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Divide(rssringoccs_ComplexDouble z0,
                                 rssringoccs_ComplexDouble z1)
{
    return z0/z1;
}
/*  End of rssringoccs_CDouble_Divide.                                  */

/*  Long double precision complex division.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Divide(rssringoccs_ComplexLongDouble z0,
                                     rssringoccs_ComplexLongDouble z1)
{
    return z0/z1;
}
/*  End of rssringoccs_CLDouble_Divide.                              */

#endif
