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
 *                           rss_ringoccs_complex                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Define complex data types and provide various functions for           *
 *      creating complex variables and performing complex arithmetic.         *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) complex.h (if your compiler supports C99. None otherwise):            *
 *      Standard library for complex types.                                   *
 ******************************************************************************
 *                                 WARNINGS                                   *
 ******************************************************************************
 *  1.) If your compiler supports C99 complex.h and you would rather use the  *
 *      built-in complex data type and complex functions with rss_ringoccs,   *
 *      you must build rss_ringoccs with such a configuration beforehand.     *
 *      See the config script config_librssringoccs.sh for details.           *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                                EXAMPLES                                    *
 ******************************************************************************
 *  Examples of all of the functions can be found in:                         *
 *      rss_ringoccs/examples/complex_examples/                               *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 13, 2020                                            *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/11/30 (Ryan Maguire):                                                *
 *      Added float and long double precision complex types.                  *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_COMPLEX_H__
#define __RSS_RINGOCCS_COMPLEX_H__

/*  Config header file containing the macros for how to build rss_ringoccs.   */
#include <rss_ringoccs/include/rss_ringoccs_config.h>

/*  Booleans defined here. Needed for the FFT and compare routines.           */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>

/*  C99 requires complex.h, but C11 makes it optional. If you have C99        *
 *  support and would like to use built-in complex numbers, change the setup  *
 *  script config_librssringoccs.sh. See that file for details.               *
 *  When tested against something simple like the complex exponential         *
 *  function, the cexp provided by glibc (GNU C library) for complex.h is     *
 *  slightly faster than rss_ringoccs. This was tested on an array of one     *
 *  million points in the complex plane, the times are as follows:            *
 *  (iMac 2017 running MacOS)                                                 *
 *      C89 Time: 0.026958                                                    *
 *      C99 Time: 0.022602                                                    *
 *  This is not to say one can't achieve better times with C89 compliant code.*
 *  rss_ringoccs uses simple, but effective algorithms, not necessary the     *
 *  fastest. Still, ~0.027 seconds for a very large array isn't bad. The code *
 *  was compiled using gcc 10 with -O3 optimization. For 10 million points    *
 *  the times were:                                                           *
 *      C89 Time: 0.267769                                                    *
 *      C99 Time: 0.231087                                                    *
 *  Which seems linear.                                                       */

/*  Check if the user has requested building using C99 compliant code.        */
#if __RSS_RINGOCCS_USING_COMPLEX_H__ == 1

/*  We'll only use C99 if the compiler supports it. The __STDC_VERSION__      *
 *  should be 199901L or higher and the __STDC_NO_COMPLEX__ macro should not  *
 *  be defined. By default the -std=c89 and -ansi flags are enabled for the   *
 *  compiler in the config_librssringoccs.sh script. This will set            *
 *  __STDC_VERSION__ to its C89 value. Remove -ansi and change -std=c89 to    *
 *  -std=c99 if you want complex.h support. If you are using a compiler other *
 *  than clang or gcc, see it's documentation for details.                    */
#if !defined(__STDC_NO_COMPLEX__) && __STDC_VERSION__ >= 199901L

/*  The various .c files will check if this is 1 or 0 to compile the right    *
 *  code, depending on if you have complex.h or not.                          */
#define _RSS_RINGOCCS_USING_COMPLEX_H_ 1

/*  Grab everything from the C99 standard complex.h and we'll just alias the  *
 *  functions and macros defined within.                                      */
#include <complex.h>

/*  You have complex.h support, so we'll just typedef double _Complex.        */
typedef double _Complex rssringoccs_ComplexDouble;

/*  Typedef single and long double precision equivalents.                     */
typedef float _Complex rssringoccs_ComplexFloat;
typedef long double _Complex rssringoccs_ComplexLongDouble;

#else
/*  Else for #if !defined(__STDC_NO_COMPLEX__) && __STDC_VERSION__ >= 199901L */

/*  If we get here, then you've requested C99 complex.h but your compiler     *
 *  not support it. Abort compiling with an error.                            */

#error __RSS_RINGOCCS_USING_COMPLEX_H__ set to 1 but you do not have complex.h

#endif
/*  End of #if !defined(__STDC_NO_COMPLEX__) && __STDC_VERSION__ >= 199901L.  */

#else
/*  Else statement for #if __RSS_RINGOCCS_USING_COMPLEX_H__ == 1.             */

/*  This is the default build for rss_ringoccs. The code is strictly C89      *
 *  compliant and attempts to be as portable as possible.                     */

/*  Set the _RSS_RINGOCCS_USING_COMPLEX_H_ macro to zero.                     */
#define _RSS_RINGOCCS_USING_COMPLEX_H_ 0

/*  The GNU Scientific Library (GSL) v2.6 defines complex variables via a     *
 *  data structure containing a single array double dat[2];. If you are using *
 *  the GSL v2.6, you can use rss_ringoccs functions with that library. That  *
 *  is, if we have a pointer rssringoccs_ComplexDouble *z; and another        *
 *  pointer gsl_complex *w; we can safely cast via:                           *
 *      z = (rssringoccs_ComplexDouble *)&w;                                  *
 *  And similarly we can do w = (gsl_complex *)&z;                            */
typedef struct rssringoccs_ComplexDouble {
    double dat[2];
} rssringoccs_ComplexDouble;

/*  Define single and long double precision equivalents.                      */
typedef struct rssringoccs_ComplexFloat {
    float dat[2];
} rssringoccs_ComplexFloat;

typedef struct rssringoccs_ComplexLongDouble {
    long double dat[2];
} rssringoccs_ComplexLongDouble;

#endif
/*  End of __RSS_RINGOCCS_USING_COMPLEX_H__ == 1.                             */

/*  Useful constants used throughout computations.                            */

/*  The imaginary unit, z = 0 + i 1.                                          */
extern const rssringoccs_ComplexFloat rssringoccs_CFloat_I;
extern const rssringoccs_ComplexDouble rssringoccs_CDouble_I;
extern const rssringoccs_ComplexLongDouble rssringoccs_CLDouble_I;

/*  Complex zero, z = 0 + i0.                                                 */
extern const rssringoccs_ComplexFloat rssringoccs_CFloat_Zero;
extern const rssringoccs_ComplexDouble rssringoccs_CDouble_Zero;
extern const rssringoccs_ComplexLongDouble rssringoccs_CLDouble_Zero;

/*  Complex one, z = 1 + i0.                                                  */
extern const rssringoccs_ComplexFloat rssringoccs_CFloat_One;
extern const rssringoccs_ComplexDouble rssringoccs_CDouble_One;
extern const rssringoccs_ComplexLongDouble rssringoccs_CLDouble_One;

/*  Complex Not-a-Number, set to NaN + i NaN.                                 */
extern const rssringoccs_ComplexFloat rssringoccs_CFloat_NaN;
extern const rssringoccs_ComplexDouble rssringoccs_CDouble_NaN;
extern const rssringoccs_ComplexLongDouble rssringoccs_CLDouble_NaN;

/*  Complex infinity, set to inf + i inf. This is the "north pole" on the     *
 *  Riemann sphere.                                                           */
extern const rssringoccs_ComplexFloat rssringoccs_CFloat_Infinity;
extern const rssringoccs_ComplexDouble rssringoccs_CDouble_Infinity;
extern const rssringoccs_ComplexLongDouble rssringoccs_CLDouble_Infinity;

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Abs                                               *
 *  Purpose:                                                                  *
 *      Compute the absolute value of a complex number. This is equivalent to *
 *      the cabs function found in complex.h (C99).                           *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      double abs_z:                                                         *
 *          The absolute value of z, computed by the Pythagorean formula. If  *
 *          z = x + iy, then abs_z = sqrt(x^2 + y^2)                          *
 ******************************************************************************/
extern float
rssringoccs_CFloat_Abs(rssringoccs_ComplexFloat z);

extern double
rssringoccs_CDouble_Abs(rssringoccs_ComplexDouble z);

extern long double
rssringoccs_CLDouble_Abs(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Abs_Squared                                       *
 *  Purpose:                                                                  *
 *      Compute the square of the absolute value of a complex number z. This  *
 *      is useful for when |z|^2 is needed so we can avoid an expensive and   *
 *      redundant square root calculation. We can just compute x*x + y*y for  *
 *      the complex number z = x + iy.                                        *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      double abs_sq:                                                        *
 *          The square of the absolute value of z, |z|^2.                     *
 ******************************************************************************/
extern float
rssringoccs_CFloat_Abs_Squared(rssringoccs_ComplexFloat z);

extern double
rssringoccs_CDouble_Abs_Squared(rssringoccs_ComplexDouble z);

extern long double
rssringoccs_CLDouble_Abs_Squared(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Add                                               *
 *  Purpose:                                                                  *
 *      Add two complex numbers.                                              *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z0:                                         *
 *          A complex number.                                                 *
 *      rssringoccs_ComplexDouble z1:                                         *
 *          Another complex number.                                           *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble sum:                                        *
 *          The sum of z0 and z1.                                             *
 *  NOTES:                                                                    *
 *      In C99, since _Complex is a built-in data type, given double _Complex *
 *      z0 and double _Complex z1, you can just do z0 + z1. In C89 we use     *
 *      structs to define complex numbers. Structs cannot be added, so we     *
 *      need a function for computing the sum of two complex values.          *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Add(rssringoccs_ComplexFloat z1,
                             rssringoccs_ComplexFloat z2);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Add(rssringoccs_ComplexDouble z1,
                              rssringoccs_ComplexDouble z2);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add(rssringoccs_ComplexLongDouble z1,
                                  rssringoccs_ComplexLongDouble z2);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Add_Imag                                          *
 *  Purpose:                                                                  *
 *      Add an imaginary number to a complex one.                             *
 *  Arguments:                                                                *
 *      double y:                                                             *
 *          An imaginary number.                                              *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble sum:                                        *
 *          The sum of y and z.                                               *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Add_Imag(float y, rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Add_Imag(double y, rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add_Imag(long double y, rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Add_Real                                          *
 *  Purpose:                                                                  *
 *      Add a real number to a complex one.                                   *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          A real number.                                                    *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble sum:                                        *
 *          The sum of x and z.                                               *
 *  NOTES:                                                                    *
 *      This function is provided for convenience. It is somewhat laborious   *
 *      to convert a real number to a complex number and then use             *
 *      rssringoccs_CDouble_Add to add the two complex numbers, so this       *
 *      function can be used to skip the intermediate step.                   *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Add_Real(float x, rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Add_Real(double x, rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add_Real(long double x, rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Argument                                          *
 *  Purpose:                                                                  *
 *      Compute the argument (phase) of a non-zero complex number. This is    *
 *      equivalent of carg found in complex.h (C99).                          *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      double arg:                                                           *
 *          The argument of z. This is the angle z makes with the positive x  *
 *          axis and is a value between -pi and pi.                           *
 *  NOTE:                                                                     *
 *      Because this function returns a value between -pi and pi, use of this *
 *      function in the square root routine returns a branch cut along the    *
 *      negative x axis.                                                      *
 *                                                                            *
 *      Using the function on the complex zero (0, 0) returns 0.0 on          *
 *      implementations that support IEEE floating-point arithmetic. This     *
 *      includes GNU's glibc/gcc and clang.                                   *
 ******************************************************************************/
extern float
rssringoccs_CFloat_Argument(rssringoccs_ComplexFloat z);

extern double
rssringoccs_CDouble_Argument(rssringoccs_ComplexDouble z);

extern long double
rssringoccs_CLDouble_Argument(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Compare                                           *
 *  Purpose:                                                                  *
 *      Compare two complex numbers z0 and z1. This returns true if both the  *
 *      real and imaginary parts of z0 and z1 are identical, and false        *
 *      otherwise.                                                            *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z0:                                         *
 *          A complex number.                                                 *
 *      rssringoccs_ComplexDouble z1:                                         *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_Bool comp:                                                *
 *          A Boolean indicating whether or not z0 and z1 are the same.       *
 ******************************************************************************/
extern rssringoccs_Bool
rssringoccs_CFloat_Compare(rssringoccs_ComplexFloat z,
                                 rssringoccs_ComplexFloat w);

extern rssringoccs_Bool
rssringoccs_CDouble_Compare(rssringoccs_ComplexDouble z,
                                  rssringoccs_ComplexDouble w);

extern rssringoccs_Bool
rssringoccs_CLDouble_Compare(rssringoccs_ComplexLongDouble z,
                                      rssringoccs_ComplexLongDouble w);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Conjugate                                         *
 *  Purpose:                                                                  *
 *      Returns the complex conjugate of a complex number z. This is          *
 *      equivalent to conj found in complex.h (C99). If z = x + iy, the       *
 *      complex conjugate of z is conj_z = x - iy.                            *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble conj_z:                                     *
 *          The complex conjugate of z.                                       *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Conjugate(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Conjugate(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Conjugate(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Cos                                               *
 *  Purpose:                                                                  *
 *      Compute the cosine of a complex number z.                             *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble cos_z:                                      *
 *          The cosine of z.                                                  *
 *  NOTE:                                                                     *
 *      We simply use the fact that cos(x+iy) = cos(x)cos(iy)-sin(x)sin(iy)   *
 *      and then invoke the definition of hyperbolic cosine and hyperbolic    *
 *      sine yielding cos(x+iy) = cos(x)cosh(y) - i * sin(x)sinh(y).          *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Cos(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Cos(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Cos(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Dist                                              *
 *  Purpose:                                                                  *
 *      Compute the distance between two points in the complex plane.         *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble cos_z:                                      *
 *          The cosine of z.                                                  *
 *  NOTE:                                                                     *
 *      This function is provided to save one from performing                 *
 *      rssringoccs_CDouble_Subtract followed by rssringoccs_CDouble_Abs.     *
 *      This saves a redundant function call and makes code look cleaner.     *
 ******************************************************************************/
extern float
rssringoccs_CFloat_Dist(rssringoccs_ComplexFloat z0,
                        rssringoccs_ComplexFloat z1);

extern double
rssringoccs_CDouble_Dist(rssringoccs_ComplexDouble z0,
                         rssringoccs_ComplexDouble z1);

extern long double
rssringoccs_CLDouble_Dist(rssringoccs_ComplexLongDouble z0,
                          rssringoccs_ComplexLongDouble z1);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Divide                                      *
 *  Purpose:                                                                  *
 *     Compute the quotient of a complex number z0 by z1.                     *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z0:                                         *
 *          A complex number.                                                 *
 *      rssringoccs_ComplexDouble z1:                                         *
 *          Another complex number.                                           *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble quotient:                                   *
 *          The complex number z0 / z1.                                       *
 *  NOTE:                                                                     *
 *      No error check is performed on whether or not z1 = 0+0i. If this is   *
 *      true, depending on your system, you will either get +infinity for both*
 *      real and imaginary parts, or an error will occur. On MacOS and Linux  *
 *      the result is NaN + iNaN.                                             *
 *                                                                            *
 *      Division is not commutative, so given (z0, z1), this returns z0/z1 and*
 *      not z1/z0. That is, we divide the first entry by the second.          *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Divide(rssringoccs_ComplexFloat z1,
                          rssringoccs_ComplexFloat z2);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Divide(rssringoccs_ComplexDouble z1,
                           rssringoccs_ComplexDouble z2);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Divide(rssringoccs_ComplexLongDouble z1,
                            rssringoccs_ComplexLongDouble z2);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Erf                                               *
 *  Purpose:                                                                  *
 *     Compute the complex error function of z, erf(z).                       *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble erf_z:                                      *
 *          The error function evaluated at z.                                *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Erf(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Erf(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Erf(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Erfc                                              *
 *  Purpose:                                                                  *
 *     Compute the complementary complex error function of z, erf(z).         *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble erf_z:                                      *
 *          The complementary error function evaluated at z.                  *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Erfc(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Erfc(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Erfc(rssringoccs_ComplexLongDouble z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Faddeeva(rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Real_Part                                         *
 *  Purpose:                                                                  *
 *      Return the real part of a complex number. This is equivalent to creal *
 *      found in complex.h (C99).                                             *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      double real:                                                          *
 *          The real part of z.                                               *
 ******************************************************************************/
extern float
rssringoccs_CFloat_Real_Part(rssringoccs_ComplexFloat z);

extern double
rssringoccs_CDouble_Real_Part(rssringoccs_ComplexDouble z);

extern long double
rssringoccs_CLDouble_Real_Part(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Imag_Part                                         *
 *  Purpose:                                                                  *
 *      Return the imaginary part of a complex number. This is equivalent to  *
 *      cimag found in complex.h (C99).                                       *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      double real:                                                          *
 *          The imaginary part of z.                                          *
 ******************************************************************************/
extern float
rssringoccs_CFloat_Imag_Part(rssringoccs_ComplexFloat z);

extern double
rssringoccs_CDouble_Imag_Part(rssringoccs_ComplexDouble z);

extern long double
rssringoccs_CLDouble_Imag_Part(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Exp                                               *
 *  Purpose:                                                                  *
 *      Compute the complex exponential of a complex number z. This is        *
 *      equivalent to cexp defined in complex.h (C99). The complex            *
 *      exponential has the same definition as the real exponential, a power  *
 *      series with terms z^n / n!.                                           *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble exp_z:                                      *
 *          The complex exponential of z.                                     *
 *  NOTE:                                                                     *
 *      The algorithm does not actually use the power series directly, but    *
 *      rather invokes Euler's formula exp(iy) = cos(y)+isin(y). Given a      *
 *      complex number z = x+iy, we have:                                     *
 *          exp(z) = exp(x + iy)                                              *
 *                 = exp(x)exp(iy)                                            *
 *                 = exp(x)cos(y) + i exp(x)sin(y)                            *
 *      So we compute using the trig functions and the real exponential.      *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Exp(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Exp(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Exp(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Sqrt                                              *
 *  Purpose:                                                                  *
 *      Compute the principal square root of a complex number. This is        *
 *      equivalent to csqrt defined in complex.h (C99).                       *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble sqrt_z:                                     *
 *          The square root of z.                                             *
 *  NOTE:                                                                     *
 *      The algorithm computes the complex square root by putting z into polar*
 *      form, z = r exp(i theta). It then returns sqrt(r) exp(i theta/2).     *
 *      This is well defined since r is non-negative. To compute theta we use *
 *      the rssringoccs_Complex_Argument function, which returns a value in   *
 *      the range -pi < theta <= pi. Because of this there is a branch cut    *
 *      along the negative x axis. rss_ringoccs does not provide the option   *
 *      to choose a different branch.                                         *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Sqrt(rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Log                                               *
 *  Purpose:                                                                  *
 *      Compute the principal log of a complex number z. This is equivalent   *
 *      to clog defined in complex.h (C99).                                   *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble ln_z:                                       *
 *          The log of z.                                                     *
 *  NOTE:                                                                     *
 *      The algorithm computes the complex log by putting z into polar        *
 *      form, z = r exp(i theta). It then returns ln(r) + i theta, where      *
 *      ln is the real valued natural log. Because of this there is a branch  *
 *      cut along the negative x axis.                                        *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, ln_z;                            *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              ln_z = rssringoccs_Complex_Log(z);                            *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(ln_z);                   *
 *              imag = rssringoccs_Complex_Imag_Part(ln_z);                   *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 0.346574 + 0.785398i. i.e. ln(1+i) = 0.3465+0.7853i. *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Log(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Log(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Log(rssringoccs_ComplexLongDouble z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Pow(rssringoccs_ComplexDouble z0,
                        rssringoccs_ComplexDouble z1);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Real_Pow(rssringoccs_ComplexDouble z, double x);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Sin                                               *
 *  Purpose:                                                                  *
 *      Compute the sine of a complex number z.                               *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble sin_z:                                      *
 *          The sine of z.                                                    *
 *  NOTE:                                                                     *
 *      We simply use the fact that sin(x+iy) = sin(x)cos(iy)+cos(x)sin(iy)   *
 *      and then invoke the definition of hyperbolic cosine and hyperbolic    *
 *      sine yielding sin(x+iy) = sin(x)cosh(y) + i * cos(x)sinh(y).          *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Sin(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Sin(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Sin(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Tan                                               *
 *  Purpose:                                                                  *
 *      Compute the tangent of a complex number z.                            *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble sin_z:                                      *
 *          The tangent of z.                                                 *
 *  NOTE:                                                                     *
 *      We compute this via tan(z) = sin(z)/cos(z) using the complex versions *
 *      of cosine and sine.                                                   *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, tan_z;                           *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              cos_z = rssringoccs_Complex_Tan(z);                           *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(tan_z);                  *
 *              imag = rssringoccs_Complex_Imag_Part(tan_z);                  *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 1.298458 + 0.634964i i.e. sin(1+i)=1.29845+0.63496i. *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Tan(rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Rect                                              *
 *  Purpose:                                                                  *
 *      Create a complex number given it's components in Cartesian format,    *
 *      also known as rectangular format. That is, given (x, y), return x+iy. *
 *  NOTE:                                                                     *
 *      In C99 you can simply do double _Complex z = x + _Complex_I*y since   *
 *      complex variables are primitive data types, but in C89 we need to     *
 *      create a struct for them (as above). Structs can't be added, so we    *
 *      need a function for creating a complex number from two doubles.       *
 *                                                                            *
 *      By default, rss_ringoccs assumes you do NOT have a C99 compliant      *
 *      compiler, hence the need for this function.                           *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          The real component of a complex number z.                         *
 *      double y:                                                             *
 *          The imaginary component of a complex number z.                    *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble z:                                          *
 *          The complex number x + iy.                                        *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Rect(float x, float y);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Rect(double x, double y);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Rect(long double x, long double y);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Polar                                             *
 *  Purpose:                                                                  *
 *      Create a complex number given it's components in polar format. That   *
 *      is, given (r, theta), return r*exp(i * theta).                        *
 *  Arguments:                                                                *
 *      double r:                                                             *
 *          A non-negative real number. The magnitude of z.                   *
 *      double theta:                                                         *
 *          The argument of z.                                                *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble z:                                          *
 *          The complex number r exp(i theta).                                *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Polar(double r, double theta);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Subtract                                          *
 *  Purpose:                                                                  *
 *      Subtract two complex numbers.                                         *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z0:                                         *
 *          A complex number.                                                 *
 *      rssringoccs_ComplexDouble z1:                                         *
 *          Another complex number.                                           *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble diff:                                       *
 *          The difference of z0 and z1, z0 - z1.                             *
 *  NOTE:                                                                     *
 *      Subtraction is not commutative, so given (z0, z1), this computes      *
 *      the first entry minus the second. That is, z0 - z1.                   *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract(rssringoccs_ComplexFloat z1,
                                  rssringoccs_ComplexFloat z2);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract(rssringoccs_ComplexDouble z1,
                                   rssringoccs_ComplexDouble z2);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract(rssringoccs_ComplexLongDouble z1,
                                       rssringoccs_ComplexLongDouble z2);

extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract_Real(float x, rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract_Real(double x, rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract_Real(long double x,
                                            rssringoccs_ComplexLongDouble z);

extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Subtract_Imag(float y, rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Subtract_Imag(double y, rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Subtract_Imag(long double y,
                                            rssringoccs_ComplexLongDouble z);



/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Multiply                                          *
 *  Purpose:                                                                  *
 *      Mutliply two complex numbers.                                         *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z0:                                         *
 *          A complex number.                                                 *
 *      rssringoccs_ComplexDouble z1:                                         *
 *          Another complex number.                                           *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble prod:                                       *
 *          The product of z0 and z1.                                         *
 *  NOTE:                                                                     *
 *      In C99, since _Complex is a built-in data type, given double _Complex *
 *      z1 and double _Complex z2, you can just do z1 * z2. In C89 we use     *
 *      structs to define complex numbers. Structs cannot be multiplied, so   *
 *      we need a function for computing the product of two complex values.   *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Multiply(rssringoccs_ComplexFloat z1,
                                  rssringoccs_ComplexFloat z2);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Multiply(rssringoccs_ComplexDouble z1,
                                   rssringoccs_ComplexDouble z2);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Multiply(rssringoccs_ComplexLongDouble z1,
                                       rssringoccs_ComplexLongDouble z2);

extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Multiply_Real(float x, rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Multiply_Real(double x, rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Multiply_Real(long double x,
                                            rssringoccs_ComplexLongDouble z);

extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Multiply_Imag(float x, rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Multiply_Imag(double x, rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Multiply_Imag(long double x,
                                            rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_CDouble_Reciprocal                                  *
 *  Purpose:                                                                  *
 *     Compute the reciprocal (or inverse) of a complex number.               *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble rcpr_z:                                     *
 *          The reciprocal of z.                                              *
 *  NOTE:                                                                     *
 *      No error check is performed on whether or not z is 0+0i. If this is   *
 *      true, depending on your system, you will either get +infinity for both*
 *      real and imaginary parts, or an error will occur. On MacOS and Linux  *
 *      the result is NaN + iNaN.                                             *
 ******************************************************************************/
extern rssringoccs_ComplexFloat
rssringoccs_CFloat_Reciprocal(rssringoccs_ComplexFloat z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Reciprocal(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Reciprocal(rssringoccs_ComplexLongDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Poly_Real_Coeffs                                  *
 *  Purpose:                                                                  *
 *      Given a set of "degree+1" number of real coefficients and a complex   *
 *      number z, computes the polynomial f(z) = a_0 + a_1 z + ... + a_N z^N  *
 *      where N = degree.                                                     *
 *  Arguments:                                                                *
 *      double *coeffs:                                                       *
 *          A pointer to an array of doubles, the polynomial coefficients.    *
 *      unsigned int degree:                                                  *
 *          The degree of the polynomial.                                     *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  NOTE:                                                                     *
 *      One error check is performed to ensure the input coeffs pointer isn't *
 *      NULL. An exit(0) will occur if it is, crashing the program. This is   *
 *      to avoid the cryptic "segmentation fault" message that occurs when    *
 *      one tries to access a NULL pointer.                                   *
 *                                                                            *
 *      It is the users responsibility to ensure coeffs points to at least    *
 *      degree + 1 double values. Since a polynomial of degree N has N+1      *
 *      coefficients (since we need to count the zeroth term a_0), coeffs     *
 *      must be of size degree+1 or more. If not, a segmentation fault may    *
 *      occur as the result of trying to access memory we don't have.         *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble poly:                                       *
 *          The complex number given by the polynomial                        *
 *          a_0 + a_1 z + ... a_N z^N                                         *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double coeffs[31];                                            *
 *              double factorial, real, imag;                                 *
 *              unsigned int n;                                               *
 *              unsigned long N = 30;                                         *
 *              rssringoccs_ComplexDouble z, poly;                            *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 0.0);                       *
 *              factorial = 1.0;                                              *
 *                                                                            *
 *              for (n=0; n<=N; ++n)                                          *
 *              {                                                             *
 *                  coeffs[n] = 1/factorial;                                  *
 *                  factorial *= n+1;                                         *
 *              }                                                             *
 *                                                                            *
 *              poly = rssringoccs_Complex_Poly_Real_Coeffs(coeffs, N, z);    *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(poly);                   *
 *              imag = rssringoccs_Complex_Imag_Part(poly);                   *
 *                                                                            *
 *              printf("f(z) = %f + i%f\n", real, imag);                      *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This output is f(z) = 2.718282 + 0.0i. Note that this polynomial is   *
 *      is just the 30th order Taylor approximation of exp(z), so the output  *
 *      is roughly exp(1) = e.                                                *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Poly_Real_Coeffs(double *coeffs, unsigned int degree,
                                           rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Poly_Complex_Coeffs(rssringoccs_ComplexDouble *coeffs,
                                              unsigned int degree,
                                              rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexDouble
rssringoccs_CDouble_Poly_Deriv_Real_Coeffs(double *coeffs,
                                                 unsigned int degree,
                                                 unsigned int deriv,
                                                 rssringoccs_ComplexDouble z);

#endif
/*  End of include guard.                                                     */
