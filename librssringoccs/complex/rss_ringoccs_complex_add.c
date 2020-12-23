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
 *                        rss_ringoccs_complex_add                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex addition.                        *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CFloat_Add:                                               *
 *      rssringoccs_CDouble_Add:                                              *
 *      rssringoccs_CLDouble_Add:                                             *
 *  Purpose:                                                                  *
 *      Adds two complex numbers:                                             *
 *                                                                            *
 *          z + w = (a + ib) + (c + id) = (a + c) + i(b + d)                  *
 *                                                                            *
 *  Arguments:                                                                *
 *      z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          A complex number.                                                 *
 *      w (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          Another complex number.                                           *
 *  Output:                                                                   *
 *      sum (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):       *
 *          The sum of z and w.                                               *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code is compatible with the C89/C90 standard. The setup script that  *
 *  is used to compile this in config_librssringoccs.sh uses gcc and has the  *
 *  -pedantic and -std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 30, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/02 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  In C99, since _Complex is a built-in data type, given double _Complex z1  *
 *  and double _Complex z2, you can just do z1 + z2. With C89 we use structs  *
 *  to define complex numbers. Structs cannot be added, so we need a function *
 *  for computing the sum of two complex values.                              */

/*  Single precision complex addition.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Add(rssringoccs_ComplexFloat z0, rssringoccs_ComplexFloat z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat sum;
    float real0, real1;
    float imag0, imag1;
    float sum_re, sum_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_CFloat_Real_Part(z0);
    real1 = rssringoccs_CFloat_Real_Part(z1);
    imag0 = rssringoccs_CFloat_Imag_Part(z0);
    imag1 = rssringoccs_CFloat_Imag_Part(z1);

    /*  The sum of two complex numbers simply adds their components.          */
    sum_re = real0 + real1;
    sum_im = imag0 + imag1;

    /*  Create the output from sum_re and sum_im and return.                  */
    sum = rssringoccs_CFloat_Rect(sum_re, sum_im);
    return sum;
}
/*  End of rssringoccs_CFloat_Add.                                            */

/*  Double precision complex addition.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Add(rssringoccs_ComplexDouble z0,
                        rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sum;
    double real0, real1;
    double imag0, imag1;
    double sum_re, sum_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_CDouble_Real_Part(z0);
    real1 = rssringoccs_CDouble_Real_Part(z1);
    imag0 = rssringoccs_CDouble_Imag_Part(z0);
    imag1 = rssringoccs_CDouble_Imag_Part(z1);

    /*  The sum of two complex numbers simply adds their components.          */
    sum_re = real0 + real1;
    sum_im = imag0 + imag1;

    /*  Create the output from sum_re and sum_im and return.                  */
    sum = rssringoccs_CDouble_Rect(sum_re, sum_im);
    return sum;
}
/*  End of rssringoccs_CDouble_Add.                                           */

/*  Long double precision complex addition.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add(rssringoccs_ComplexLongDouble z0,
                         rssringoccs_ComplexLongDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble sum;
    long double real0, real1;
    long double imag0, imag1;
    long double sum_re, sum_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_CLDouble_Real_Part(z0);
    real1 = rssringoccs_CLDouble_Real_Part(z1);
    imag0 = rssringoccs_CLDouble_Imag_Part(z0);
    imag1 = rssringoccs_CLDouble_Imag_Part(z1);

    /*  The sum of two complex numbers simply adds their components.          */
    sum_re = real0 + real1;
    sum_im = imag0 + imag1;

    /*  Create the output from sum_re and sum_im and return.                  */
    sum = rssringoccs_CLDouble_Rect(sum_re, sum_im);
    return sum;
}
/*  End of rssringoccs_CLDouble_Add.                                          */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here we have complex.h support so we'll just use the            *
 *  arithmetic made available in the C99 standard.                            */

/*  Single precision complex addition.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Add(rssringoccs_ComplexFloat z0, rssringoccs_ComplexFloat z1)
{
    return z0 + z1;
}
/*  End of rssringoccs_CFloat_Add.                                            */

/*  Double precision complex addition.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Add(rssringoccs_ComplexDouble z0,
                        rssringoccs_ComplexDouble z1)
{
    return z0 + z1;
}
/*  End of rssringoccs_CDouble_Add.                                           */

/*  Long double precision complex addition.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add(rssringoccs_ComplexLongDouble z0,
                         rssringoccs_ComplexLongDouble z1)
{
    return z0 + z1;
}
/*  End of rssringoccs_CLDouble_Add.                                          */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */
