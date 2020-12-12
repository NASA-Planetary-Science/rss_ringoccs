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
 *                     rss_ringoccs_complex_add_real                          *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex addition.                        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 3, 2020                                              *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/03 (Ryan Maguire):                                                *
 *      Moved here from rss_ringoccs_complex_add.c.                           *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  In C99, since _Complex is a built-in data type, doubles and _Complex      *
 *  doubles can be added via x + z. With C89 we use structs to define complex *
 *  numbers. Since we can't add a double to a struct, we need a function      *
 *  for computing the sum of complex numbers with real ones.                  */

/*  Single precision complex addition where one variable is real.             */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Add_Real(float x, rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat sum;
    float real, imag, sum_re;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CFloat_Real_Part(z);
    imag = rssringoccs_CFloat_Imag_Part(z);

    /*  Add the real variable to the real part of z.                          */
    sum_re = x + real;

    /*  Create the output from sum_re and imag and return.                    */
    sum = rssringoccs_CFloat_Rect(sum_re, imag);
    return sum;
}
/*  End of rssringoccs_CFloat_Add_Real.                                       */

/*  Double precision complex addition where one variable is real.             */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Add_Real(double x, rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sum;
    double real, imag, sum_re;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CDouble_Real_Part(z);
    imag = rssringoccs_CDouble_Imag_Part(z);

    /*  Add the real variable to the real part of z.                          */
    sum_re = x + real;

    /*  Create the output from sum_re and imag and return.                    */
    sum = rssringoccs_CDouble_Rect(sum_re, imag);
    return sum;
}
/*  End of rssringoccs_CDouble_Add_Real.                                      */

/*  Long double precision complex addition where one variable is real.        */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add_Real(long double x, rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble sum;
    long double real, imag, sum_re;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CLDouble_Real_Part(z);
    imag = rssringoccs_CLDouble_Imag_Part(z);

    /*  Add the real variable to the real part of z.                          */
    sum_re = x + real;

    /*  Create the output from sum_re and imag and return.                    */
    sum = rssringoccs_CLDouble_Rect(sum_re, imag);
    return sum;
}
/*  End of rssringoccs_CLDouble_Add_Real.                                     */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here we have complex.h support so we'll just use the            *
 *  arithmetic made available in the C99 standard.                            */

/*  Single precision complex addition.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Add_Real(float x, rssringoccs_ComplexFloat z)
{
    return x + z;
}
/*  End of rssringoccs_CFloat_Add_Real.                                       */

/*  Double precision complex addition.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Add_Real(double x, rssringoccs_ComplexDouble z)
{
    return x + z;
}
/*  End of rssringoccs_CDouble_Add_Real.                                      */

/*  Long double precision complex addition.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add_Real(long double x, rssringoccs_ComplexLongDouble z)
{
    return x + z;
}
/*  End of rssringoccs_CLDouble_Add_Real.                                     */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */