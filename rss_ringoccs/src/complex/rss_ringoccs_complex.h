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
 *      Define rssringoccs_ComplexDouble data type and provide various        *
 *      functions for creating complex variables and performing basic         *
 *      arithmetic. If complex.h is available, we simply alias the functions. *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 *  1.) complex.h (if your compiler supports C99. None otherwise):            *
 *      Standard library for complex types.                                   *
 ******************************************************************************
 *                                 WARNINGS                                   *
 *  1.) If your compiler supports C99 complex.h and you would rather use the  *
 *      built-in complex data type and complex functions with rss_ringoccs,   *
 *      you must uncomment #define __RSS_RINGOCCS_HAS_COMPLEX_H_.             *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 13, 2020                                            *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_COMPLEX_H
#define RSS_RINGOCCS_COMPLEX_H

/*  Booleans defined here. Needed for the FFT routines.                       */
#include <rss_ringoccs/src/rss_ringoccs_bool.h>

/*  This is commented out since rss_ringoccs does not assume C99 compliance.  */

/*
 * #define __RSS_RINGOCCS_HAS_COMPLEX_H_
 */

/*  C99 requires complex.h, but C11 makes it optional. If you have C99        *
 *  support and would like to use built-in complex numbers, uncomment out     *
 *  #define __RSS_RINGOCCS_HAS_COMPLEX_H_ above, then rebuild rss_ringoccs.   *
 *  When tested against something simple like the complex exponential         *
 *  function, the cexp provided by glibc (GNU C library) for complex.h is     *
 *  slightly faster than rss_ringoccs. This was tested on an array of one     *
 *  million points in the complex plane, the times are as follows:            *
 *  (iMac 2017 running MacOS)                                                 *
 *      C89 Time: 0.026958                                                    *
 *      C99 Time: 0.022602                                                    *
 *  This is not to say one can't achieve better times with C89 compliant code.*
 *  rss_ringoccs uses simple, but effective algorithms, not necessary the     *
 *  fastest. Still, ~0.027 seconds for a very large array isn't bad.          *
 *  Note, the two times computed above use rss_ringoccs with C89 support and  *
 *  rss_ringoccs with C99 support, respectively. rss_ringoccs simply aliases  *
 *  cexp as rssringoccs_Complex_Exp, so they compute the exact same thing,    *
 *  whereas C89 rss_ringoccs has an algorithm implemented for the complex     *
 *  exponential. The code was compiled using gcc 10 with -O3 optimization.    *
 *  For 10 million points the times were:                                     *
 *      C89 Time: 0.267769                                                    *
 *      C99 Time: 0.231087                                                    *
 *  Which seems linear.                                                       */

/*  We'll only use C99 if the user requested (defined the                     *
 *  __RSS_RINGOCCS_HAS_COMPLEX_H_ macro), and the compiler supports it. The   *
 *  __STDC_VERSION__ should be, at least, C99 capable, and the                *
 *  __STDC_NO_COMPLEX__ macro should not be defined.                          */
#if defined(__RSS_RINGOCCS_HAS_COMPLEX_H_) &&                                  \
    !defined(__STDC_NO_COMPLEX__)          &&                                  \
    __STDC_VERSION__ >= 199901L

/*  The various .c files will check if this is 1 or 0 to compile the right    *
 *  code, depending on if you have complex.h or not.                          */
#define _RSS_RINGOCCS_USING_COMPLEX_H_ 1

/*  Grab everything from the C99 standard complex.h and we'll just alias the  *
 *  functions and macros defined within.                                      */
#include <complex.h>

/*  You have complex.h support, so we'll just alias the various functions.    */
typedef double _Complex rssringoccs_ComplexDouble;
#define rssringoccs_Complex_Abs cabs
#define rssringoccs_Complex_Argument carg
#define rssringoccs_Complex_Real_Part creal
#define rssringoccs_Complex_Imag_Part cimag
#define rssringoccs_Complex_Conjugate conj
#define rssringoccs_Complex_Exp cexp
#define rssringoccs_Complex_Sqrt csqrt
#define rssringoccs_Complex_Log clog
#define rssringoccs_Complex_Sin csin
#define rssringoccs_Complex_Cos ccos
#define rssringoccs_Complex_Tan ctan

#else
/*  If we get here, your compiler does not support complex.h, or you have     *
 *  chosen to use C89 compliant code. Either way, we'll need to create        *
 *  various functions to enable complex arithmetic.                           */

/*  Set the _RSS_RINGOCCS_USING_COMPLEX_H_ macro to zero.                     */
#define _RSS_RINGOCCS_USING_COMPLEX_H_ 0

/*  The GNU Scientific Library (GSL) v2.6 defines complex variables via a     *
 *  data structure containing a single array double dat[2];. If you are using *
 *  the GSL v2.6, you can use rss_ringoccs functions with that library.       */
typedef struct {
    double dat[2];
} rssringoccs_ComplexDouble;

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Abs                                               *
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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double abs_z;                                                 *
 *              rssringoccs_ComplexDouble z;                                  *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              abs_z = rssringoccs_Complex_Abs(z);                           *
 *                                                                            *
 *              printf("%f\n", abs_z);                                        *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 1.414..., the square root of two.                    *
 ******************************************************************************/
extern double rssringoccs_Complex_Abs(rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Argument                                          *
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
 *      included GNU's glibc/gcc and clang.                                   *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double theta;                                                 *
 *              rssringoccs_ComplexDouble z;                                  *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              theta = rssringoccs_Complex_Argument(z);                      *
 *                                                                            *
 *              printf("%f\n", theta);                                        *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 0.7853..., pi/4.                                     *
 ******************************************************************************/
extern double rssringoccs_Complex_Argument(rssringoccs_ComplexDouble z);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double z_real;                                                *
 *              rssringoccs_ComplexDouble z;                                  *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              z_real = rssringoccs_Complex_Real_Part(z);                    *
 *                                                                            *
 *              printf("%f\n", z_real);                                       *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 1.0.                                                 *
 ******************************************************************************/
extern double rssringoccs_Complex_Real_Part(rssringoccs_ComplexDouble z);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double z_imag;                                                *
 *              rssringoccs_ComplexDouble z;                                  *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              z_imag = rssringoccs_Complex_Imag_Part(z);                    *
 *                                                                            *
 *              printf("%f\n", z_imag);                                       *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 1.0.                                                 *
 ******************************************************************************/
extern double rssringoccs_Complex_Imag_Part(rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Conjugate                                         *
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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, conj_z;                          *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              conj_z = rssringoccs_Complex_Conjugate(z);                    *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(conj_z);                 *
 *              imag = rssringoccs_Complex_Imag_Part(conj_z);                 *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 1.0 + -1.0i.                                         *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Conjugate(rssringoccs_ComplexDouble z);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_math.h>                                    *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, exp_z;                           *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(0.0, ONE_PI);                    *
 *              exp_z = rssringoccs_Complex_Exp(z);                           *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(exp_z);                  *
 *              imag = rssringoccs_Complex_Imag_Part(exp_z);                  *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut -1.0 + 0.0i. i.e., exp(i pi) = -1                    *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Exp(rssringoccs_ComplexDouble z);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_math.h>                                    *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, sqrt_z;                          *
 *                                                                            *
 *              sqrt_z = rssringoccs_Complex_Sqrt(rssringoccs_Imaginary_Unit);*
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(sqrt_z);                 *
 *              imag = rssringoccs_Complex_Imag_Part(sqrt_z);                 *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(-1.0, 0.0);                      *
 *              sqrt_z = rssringoccs_Complex_Sqrt(z);                         *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(sqrt_z);                 *
 *              imag = rssringoccs_Complex_Imag_Part(sqrt_z);                 *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 0.707107 + 0.707107i, i.e. (1+i)/sqrt(2), and then   *
 *      0.000000 + 1.000000i. This exemplifies the branch cut. We have        *
 *      sqrt(-1) = 0 + 1i, and not 0 + -1i.                                   *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Sqrt(rssringoccs_ComplexDouble z);

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
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Log(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexDouble
rssringoccs_Complex_Pow(rssringoccs_ComplexDouble z0,
                        rssringoccs_ComplexDouble z1);

extern rssringoccs_ComplexDouble
rssringoccs_Complex_Real_Pow(rssringoccs_ComplexDouble z, double x);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, sin_z;                           *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              sin_z = rssringoccs_Complex_Sin(z);                           *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(sin_z);                  *
 *              imag = rssringoccs_Complex_Imag_Part(sin_z);                  *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 1.298458 + 0.634964i i.e. sin(1+i)=1.29845+0.63496i. *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Sin(rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Cos                                               *
 *  Purpose:                                                                  *
 *      Compute the cosine of a complex number z.                             *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble sin_z:                                      *
 *          The cosine of z.                                                  *
 *  NOTE:                                                                     *
 *      We simply use the fact that cos(x+iy) = cos(x)cos(iy)-sin(x)sin(iy)   *
 *      and then invoke the definition of hyperbolic cosine and hyperbolic    *
 *      sine yielding cos(x+iy) = cos(x)cosh(y) - i * sin(x)sinh(y).          *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, cos_z;                           *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              cos_z = rssringoccs_Complex_Cos(z);                           *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(cos_z);                  *
 *              imag = rssringoccs_Complex_Imag_Part(cos_z);                  *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This will ouptut 1.298458 + 0.634964i i.e. sin(1+i)=1.29845+0.63496i. *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Cos(rssringoccs_ComplexDouble z);

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

#endif
/*  End of #ifdef __RSS_RINGOCCS_HAS_COMPLEX_H_                               */

/*  Useful constants used throughout computations.                            */
extern const rssringoccs_ComplexDouble rssringoccs_Imaginary_Unit;
extern const rssringoccs_ComplexDouble rssringoccs_Complex_Zero;
extern const rssringoccs_ComplexDouble rssringoccs_Complex_One;

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double abs_z, real, imag;                                     *
 *              rssringoccs_ComplexDouble z;                                  *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              abs_z = rssringoccs_Complex_Abs(z);                           *
 *              real = rssringoccs_Complex_Real_Part(z);                      *
 *              imag = rssringoccs_Complex_Imag_Part(z);                      *
 *                                                                            *
 *              printf("|%f + %fi| = %f\n", real, imag, abs_z);               *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This creates the complex number 1+1i. The output prints:              *
 *          |1.000000 + 1.000000i| = 1.414214                                 *
 ******************************************************************************/
extern rssringoccs_ComplexDouble rssringoccs_Complex_Rect(double x, double y);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_math.h>                                    *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z;                                  *
 *                                                                            *
 *              z = rssringoccs_Complex_Polar(1.0, PI_BY_FOUR);               *
 *              real = rssringoccs_Complex_Real_Part(z);                      *
 *              imag = rssringoccs_Complex_Imag_Part(z);                      *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs 0.707107 + 0.707107i, i.e. (1+i1)/sqrt(2).               *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Polar(double r, double theta);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Add                                               *
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
 *  NOTE:                                                                     *
 *      In C99, since _Complex is a built-in data type, given double _Complex *
 *      z0 and double _Complex z1, you can just do z0 + z1. In C89 we use     *
 *      structs to define complex numbers. Structs cannot be added, so we     *
 *      need a function for computing the sum of two complex values.          *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z0, z1, sum;                        *
 *                                                                            *
 *              z0 = rssringoccs_Complex_Rect(1.0, 2.0);                      *
 *              z1 = rssringoccs_Complex_Rect(2.0, 2.0);                      *
 *                                                                            *
 *              sum = rssringoccs_Complex_Add(z0, z1);                        *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(sum);                    *
 *              imag = rssringoccs_Complex_Imag_Part(sum);                    *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs 3.0 + 4.0i, i.e. (1+2i) + (2+2i) = 3+4i.                 *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Add(rssringoccs_ComplexDouble z1,
                        rssringoccs_ComplexDouble z2);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z0, z1, diff;                       *
 *                                                                            *
 *              z0 = rssringoccs_Complex_Rect(1.0, 2.0);                      *
 *              z1 = rssringoccs_Complex_Rect(2.0, 2.0);                      *
 *                                                                            *
 *              diff = rssringoccs_Complex_Subtract(z0, z1);                  *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(diff);                   *
 *              imag = rssringoccs_Complex_Imag_Part(diff);                   *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs -1.0 + 0.0i, i.e. (1+2i) - (2+2i) = -1+0i.               *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Subtract(rssringoccs_ComplexDouble z1,
                             rssringoccs_ComplexDouble z2);

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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z0, z1, prod;                       *
 *                                                                            *
 *              z0 = rssringoccs_Complex_Rect(1.0, 2.0);                      *
 *              z1 = rssringoccs_Complex_Rect(2.0, 2.0);                      *
 *                                                                            *
 *              prod = rssringoccs_Complex_Multiply(z0, z1);                  *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(prod);                   *
 *              imag = rssringoccs_Complex_Imag_Part(prod);                   *
 *                                                                            *
 *              printf("%f + %fi\n", real, imag);                             *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs -2.0 + 6.0i, i.e. (1+2i)*(2+2i) = -2+6i.                 *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Multiply(rssringoccs_ComplexDouble z1,
                             rssringoccs_ComplexDouble z2);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Scale                                             *
 *  Purpose:                                                                  *
 *     Scale a complex number by a real one.                                  *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          A real number.                                                    *
 *      rssringoccs_ComplexDouble z1:                                         *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble scale:                                      *
 *          The product x*z.                                                  *
 *  NOTE:                                                                     *
 *      In C99, since _Complex is a built-in data type, given double x and a  *
 *      double _Complex z, you can just do x * z. In C89 we use structs to    *
 *      define complex numbers. Structs cannot be multiplied by real numbers  *
 *      so we need a function for computing the product of a real number with *
 *      a complex one.                                                        *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double x, abs_z, abs_scale;                                   *
 *              rssringoccs_ComplexDouble z, scale;                           *
 *                                                                            *
 *              x = 2.0;                                                      *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              scale = rssringoccs_Complex_Scale(x, z);                      *
 *                                                                            *
 *              abs_z = rssringoccs_Complex_Abs(z);                           *
 *              abs_scale = rssringoccs_Complex_Abs(scale);                   *
 *                                                                            *
 *              printf("|z| = %f\t |xz| = %f\n", abs_z, abs_scale);           *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs |z| = 1.414214	 |xz| = 2.828427. i.e., |z| = sqrt(2) and *
 *      |2z| = 2*sqrt(2).                                                     *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Scale(double x, rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Reciprocal                                        *
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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z, inv_z;                           *
 *                                                                            *
 *              z = rssringoccs_Complex_Rect(1.0, 1.0);                       *
 *              inv_z = rssringoccs_Complex_Reciprocal(z);                    *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(inv_z);                  *
 *              imag = rssringoccs_Complex_Imag_Part(inv_z);                  *
 *                                                                            *
 *              printf("1/z = %f + i%f\n", real, imag);                       *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs 1/z = 0.5 + i-0.5. i.e., 1/(1+1i) = 0.5 - 0.5i.          *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Reciprocal(rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Divide                                            *
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
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z0, z1, quotient;                   *
 *                                                                            *
 *              z0 = rssringoccs_Complex_Rect(2.0, 1.0);                      *
 *              z1 = rssringoccs_Complex_Rect(-1.0, 4.0);                     *
 *              quotient = rssringoccs_Complex_Divide(z0, z1);                *
 *                                                                            *
 *              real = rssringoccs_Complex_Real_Part(quotient);               *
 *              imag = rssringoccs_Complex_Imag_Part(quotient);               *
 *                                                                            *
 *              printf("z0/z1 = %f + i%f\n", real, imag);                     *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs z0/z1 = 0.117647 + i-0.529412.                           *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Divide(rssringoccs_ComplexDouble z1,
                           rssringoccs_ComplexDouble z2);

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
rssringoccs_Complex_Poly_Real_Coeffs(double *coeffs, unsigned int degree,
                                     rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexDouble
rssringoccs_Complex_Poly_Complex_Coeffs(rssringoccs_ComplexDouble *coeffs,
                                        unsigned int degree,
                                        rssringoccs_ComplexDouble z);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Complex_Erf                                               *
 *  Purpose:                                                                  *
 *     Compute the complex error function of z, erf(z).                       *
 *  Arguments:                                                                *
 *      rssringoccs_ComplexDouble z:                                          *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      rssringoccs_ComplexDouble erf_z:                                      *
 *          The error function evaluated at z.                                *
 *  NOTE:                                                                     *
 *      The algorithm implemented uses a Taylor series for complex values     *
 *      with magnitude less than 3 and an asymptotic series otherwise. The    *
 *      error is about 10^-7 for values within a reasonable distance from the *
 *      origin (say, magnitude less than 10). For certain arguments the       *
 *      error function quickly converges to 1, and for others it explodes     *
 *      asymptotically like ~exp(x^2).                                        *
 *  Example:                                                                  *
 *          #include <stdio.h>                                                *
 *          #include <rss_ringoccs_complex.h>                                 *
 *                                                                            *
 *          int main(void)                                                    *
 *          {                                                                 *
 *              double real, imag;                                            *
 *              rssringoccs_ComplexDouble z0, z1, erf_z;                      *
 *                                                                            *
 *              z0 = rssringoccs_Complex_Rect(1.0, 1.0);                      *
 *              z1 = rssringoccs_Complex_Rect(10.0, 4.0);                     *
 *                                                                            *
 *              erf_z = rssringoccs_Complex_Erf(z0);                          *
 *              real = rssringoccs_Complex_Real_Part(erf_z);                  *
 *              imag = rssringoccs_Complex_Imag_Part(erf_z);                  *
 *              printf("erf(1+1i) = %f + i%f\n", real, imag);                 *
 *                                                                            *
 *              erf_z = rssringoccs_Complex_Erf(z1);                          *
 *              real = rssringoccs_Complex_Real_Part(erf_z);                  *
 *              imag = rssringoccs_Complex_Imag_Part(erf_z);                  *
 *              printf("erf(10+4i) = %f + i%f\n", real, imag);                *
 *                                                                            *
 *              return 0;                                                     *
 *          }                                                                 *
 *                                                                            *
 *      This outputs erf(1+1i) = 1.316151 + i0.190453 and then                *
 *      erf(10+4i) = 1.000000 + i-0.000000, demonstrating how quickly erf     *
 *      converges to 1 for certain arguments (angles/phases).                 *
 ******************************************************************************/
extern rssringoccs_ComplexDouble
rssringoccs_Complex_Erf(rssringoccs_ComplexDouble z);

extern rssringoccs_ComplexDouble *
rssringoccs_FFT_Cooley_Tukey_ComplexDouble(rssringoccs_ComplexDouble *in,
                                           long N, rssringoccs_Bool inverse);

extern rssringoccs_ComplexDouble *
rssringoccs_FFT_Bluestein_Chirp_Z_ComplexDouble(rssringoccs_ComplexDouble *in,
                                                long N,
                                                rssringoccs_Bool inverse);

#endif
/*  End of include guard.                                                     */
