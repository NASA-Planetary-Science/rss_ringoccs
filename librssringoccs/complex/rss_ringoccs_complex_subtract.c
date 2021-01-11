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
 *                      rss_ringoccs_complex_subtract                         *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex subtraction.                     *
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
 *  and double _Complex z2, you can just do z1 + z2. Structs can't be added,  *
 *  so we need a function for computing the sum of two complex values.        */

/*  Single precision complex addition.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract(rssringoccs_ComplexFloat z0,
                                  rssringoccs_ComplexFloat z1)
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
    sum_re = real0 - real1;
    sum_im = imag0 - imag1;

    /*  Create the output from sum_re and sum_im and return.                  */
    sum = rssringoccs_CFloat_Rect(sum_re, sum_im);
    return sum;
}

/*  Double precision complex addition.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract(rssringoccs_ComplexDouble z0,
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
    sum_re = real0 - real1;
    sum_im = imag0 - imag1;

    /*  Create the output from sum_re and sum_im and return.                  */
    sum = rssringoccs_CDouble_Rect(sum_re, sum_im);
    return sum;
}

/*  Long double precision complex addition.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract(rssringoccs_ComplexLongDouble z0,
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
    sum_re = real0 - real1;
    sum_im = imag0 - imag1;

    /*  Create the output from sum_re and sum_im and return.                  */
    sum = rssringoccs_CLDouble_Rect(sum_re, sum_im);
    return sum;
}

/*  Single precision complex addition where one variable is real.             */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract_Real(float x, rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat sum;
    float real, imag, sum_re;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CFloat_Real_Part(z);
    imag = rssringoccs_CFloat_Imag_Part(z);

    /*  Add the real variable to the real part of z.                          */
    sum_re = x - real;

    /*  Create the output from sum_re and imag and return.                    */
    sum = rssringoccs_CFloat_Rect(sum_re, -imag);
    return sum;
}

/*  Double precision complex addition where one variable is real.             */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract_Real(double x, rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sum;
    double real, imag, sum_re;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CDouble_Real_Part(z);
    imag = rssringoccs_CDouble_Imag_Part(z);

    /*  Add the real variable to the real part of z.                          */
    sum_re = x - real;

    /*  Create the output from sum_re and imag and return.                    */
    sum = rssringoccs_CDouble_Rect(sum_re, -imag);
    return sum;
}

/*  Long double precision complex addition where one variable is real.        */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract_Real(long double x,
                                            rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble sum;
    long double real, imag, sum_re;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CLDouble_Real_Part(z);
    imag = rssringoccs_CLDouble_Imag_Part(z);

    /*  Add the real variable to the real part of z.                          */
    sum_re = x - real;

    /*  Create the output from sum_re and imag and return.                    */
    sum = rssringoccs_CLDouble_Rect(sum_re, -imag);
    return sum;
}

/*  Single precision complex addition where one variable is imaginary.        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract_Imag(float y, rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat sum;
    float real, imag, sum_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CFloat_Real_Part(z);
    imag = rssringoccs_CFloat_Imag_Part(z);

    /*  Add the imaginary variable to the imaginary part of z.                */
    sum_im = y - imag;

    /*  Create the output from real and sum_im and return.                    */
    sum = rssringoccs_CFloat_Rect(-real, sum_im);
    return sum;
}

/*  Double precision complex addition where one variable is imaginary.        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract_Imag(double y, rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sum;
    double real, imag, sum_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CDouble_Real_Part(z);
    imag = rssringoccs_CDouble_Imag_Part(z);

    /*  Add the imaginary variable to the imaginary part of z.                */
    sum_im = y - imag;

    /*  Create the output from real and sum_im and return.                    */
    sum = rssringoccs_CDouble_Rect(-real, sum_im);
    return sum;
}

/*  Long double precision complex addition where one variable is real.        */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract_Imag(long double y,
                                            rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble sum;
    long double real, imag, sum_im;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real = rssringoccs_CLDouble_Real_Part(z);
    imag = rssringoccs_CLDouble_Imag_Part(z);

    /*  Add the imaginary variable to the imaginary part of z.                */
    sum_im = y - imag;

    /*  Create the output from real and sum_im and return.                    */
    sum = rssringoccs_CLDouble_Rect(-real, sum_im);
    return sum;
}

#else
/*  If we get here we have complex.h support so we'll just alias the          *
 *  arithmetic made available with the C99 standard.                          */

/*  Single precision complex addition.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract(rssringoccs_ComplexFloat z0,
                                  rssringoccs_ComplexFloat z1)
{
    return z0 - z1;
}

rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract_Real(float x, rssringoccs_ComplexFloat z)
{
    return x - z;
}

rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract_Imag(float y, rssringoccs_ComplexFloat z)
{
    return _Complex_I*y - z;
}

/*  Double precision complex addition.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract(rssringoccs_ComplexDouble z0,
                                   rssringoccs_ComplexDouble z1)
{
    return z0 - z1;
}

rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract_Real(double x, rssringoccs_ComplexDouble z)
{
    return x - z;
}

rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract_Imag(double y, rssringoccs_ComplexDouble z)
{
    return _Complex_I*y - z;
}

/*  Long double precision complex addition.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract(rssringoccs_ComplexLongDouble z0,
                                       rssringoccs_ComplexLongDouble z1)
{
    return z0 - z1;
}

rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract_Real(long double x,
                                            rssringoccs_ComplexLongDouble z)
{
    return x - z;
}


rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract_Imag(long double y,
                                            rssringoccs_ComplexLongDouble z)
{
    return _Complex_I*y - z;
}

#endif
