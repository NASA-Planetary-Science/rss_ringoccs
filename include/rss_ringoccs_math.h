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
 *                            rss_ringoccs_math                               *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This header file is provided to improve the portability of            *
 *      rss_ringoccs. The C99 math library is a superset of the C89/C90       *
 *      version and contains float and long double support for all of the     *
 *      standard functions. So if "func" is provided in C89/C90, C99 contains *
 *      func as well as:                                                      *
 *          funcf:                                                            *
 *              Float version of func.                                        *
 *          funcl:                                                            *
 *              Long double version of func.                                  *
 *      For example, sinf, sinl, fabsf, fabsl, etc. The older C89/C90 math.h  *
 *      do not provide these functions, and rather if a function recieves a   *
 *      non-double (like passing a float to cos), then an implicit type       *
 *      conversion occurs, which may be computationally expensive. The funcf  *
 *      and funcl versions are meant to rid of this conversion step.          *
 *      For float and long double functions rss_ringoccs uses, if available,  *
 *      these functions. Here we provide aliases for the functions in math.h  *
 *      depending on whether or not the __HAS_C99_MATH_H__ macro is defined   *
 *      (it's defined in rss_ringoccs_config.h).                              *
 *                                                                            *
 *      There are also new functions in the C99 version, such as the          *
 *      hyperbolic trig functions, the copysign function, and more. We provide*
 *      the prototypes for these here, and the algorithms for these functions *
 *      can be found in the math/ subdirectory of librssringoccs.             *
 *                                                                            *
 *      We do not assume you are using a C99 capable compiler and instead     *
 *      assume strict C89/C90 (also called ANSI C) compliance. The code still *
 *      compiles on C99 and C11 compliant compilers.                          *
 *                                                                            *
 *      This file also provides NaN and infinity macros if they are not set.  *
 *      NOTE:                                                                 *
 *          Infinity is set as the standard macro HUGE_VAL defined in math.h  *
 *          and for most implementations this should do. Indeed, this is the  *
 *          same manner the Py_HUGE_VAL is set. The python source code issues *
 *          the following warning (cpython/Include/pymath.h):                 *
 *              HUGE_VAL is supposed to expand to a positive double infinity. *
 *              Python uses Py_HUGE_VAL instead because some platforms are    *
 *              broken in this respect.  We used to embed code in pyport.h to *
 *              try to worm around that, but different platforms are broken   *
 *              in conflicting ways.  If you're on a platform where HUGE_VAL  *
 *              is defined incorrectly, fiddle your Python config to          *
 *              #define Py_HUGE_VAL to something that works on your platform. *
 *                                                                            *
 *          Similarly, NaN is defined as HUGE_VAL * 0, which should be        *
 *          infinity times zero, which is Not-A-Number. Python does this as   *
 *          well, issuing the following warning:                              *
 *              Py_NAN                                                        *
 *              A value that evaluates to a NaN. On IEEE 754 platforms INF*0  *
 *              or INF/INF works. Define Py_NO_NAN in pyconfig.h if your      *
 *              platform doesn't support NaNs.                                *
 *          If necessary, redefine NAN here to whatever your platform allows. *
 *      Lastly, this file provides a bunch of constants that are commonly     *
 *      used, like sqrt(2) and multiples of pi.                               *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) math.h:                                                               *
 *      Standard library for mathematical functions like sin, cos, atan.      *
 *  2.) float.h:                                                              *
 *      Standard library which contains macros for the smallest and largest   *
 *      values allowed by your system.                                        *
 *  3.) rss_ringoccs_config.h:                                                *
 *      Header file which contains the __HAS_C99_MATH_H__ macro. This macro   *
 *      tells the compiler how to build librssringoccs.
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 12, 2020                                            *
 ******************************************************************************/

/*  Include guard for this file to prevent including it twice.                */
#ifndef __RSS_RINGOCCS_MATH_H__
#define __RSS_RINGOCCS_MATH_H__

/*  The __RSS_RINGOCCS_USING_C99_MATH_H__ macro is found here.                */
#include <rss_ringoccs/include/rss_ringoccs_config.h>

/*  It is not assumed you have C99 math.h, but rather the C89 math.h, which   *
 *  is also referred to as the ANSI C version. If you have C99 available,     *
 *  change __RSS_RINGOCCS_USING_C99_MATH_H__ in rss_ringoccs_config.h to 1    *
 *  and then rebuild rss_ringoccs. Some implementations of the C99 standard   *
 *  don't define float/long double versions of math functions and simply do   *
 *  something like sinf(x) = (float)sin((double)x), so there is no difference.*/

/*  Check which version of C you are compiling with and set the macro         *
 *  __HAS_C99_MATH_H__ accordingly.                                           */
#if __RSS_RINGOCCS_USING_C99_MATH_H__ == 1

/*  The __STDC_VERSION__ macro should be set to C99 or higher. Check this.    */
#if __STDC_VERSION__ >= 199901L

/*  You've requested C99 math.h and have support for it, so set the           *
 *  __HAS_C99_MATH_H__ macro to 1.                                            */
#define __HAS_C99_MATH_H__ 1

#else
/*  Else statement for #if __STDC_VERSION__ >= 199901L.                       */

/*  You requested C99 math.h but don't have it. Abort with error.             */
#error Requested C99 math.h but your compiler does not support it.

#endif
/*  End of #if __STDC_VERSION__ >= 199901L                                    */

#else
/*  Else statement for #if __RSS_RINGOCCS_USING_C99_MATH_H__ == 1.            */

/*  You have not requested C99 math.h, so set __HAS_C99_MATH_H__ to 0.        */
#define __HAS_C99_MATH_H__ 0

#endif
/*  End of #if __RSS_RINGOCCS_USING_C99_MATH_H__ == 1.                        */

/*  Include the standard library header math.h. We're only going to alias     *
 *  functions we ever use in rss_ringoccs, sin, cos, fabs, exp, atan2, etc.   */
#include <math.h>

/*  And this header file contains macros for the smallest and largest allowed *
 *  values for your system.                                                   */
#include <float.h>

/* Define Miscellaneous Constants.                                            */

/*  Single precision constants.                                               */
#define rssringoccs_Sqrt_One_By_Two_Pi_F    0.398942280F
#define rssringoccs_Sqrt_Pi_By_Eight_F      0.626657069F
#define rssringoccs_Sqrt_Pi_By_Two_F        1.253314137F
#define rssringoccs_Sqrt_One_By_Pi_F        0.564189584F
#define rssringoccs_Sqrt_Two_By_Pi_F        0.797884561F
#define rssringoccs_Two_By_Sqrt_Pi_F        1.128379167F
#define rssringoccs_Pi_By_Two_F             1.570796327F
#define rssringoccs_Three_Pi_By_Four_F      2.356194490F
#define rssringoccs_Pi_By_Four_F            0.785398163F
#define rssringoccs_One_Pi_F                3.141592654F
#define rssringoccs_Two_Pi_F                6.283185307F
#define rssringoccs_Sqrt_Two_F              1.414213562F
#define rssringoccs_Rcpr_Euler_E_F          0.367879441F
#define rssringoccs_Euler_E_F               2.718281828F
#define rssringoccs_Natural_Log_of_2_F      0.693147181F
#define rssringoccs_Natural_Log_Of_10_F     2.302585093F

/*  Double precision constants.                                               */
#define rssringoccs_Sqrt_One_By_Two_Pi      0.39894228040143267
#define rssringoccs_Sqrt_Pi_By_Eight        0.62665706865775012
#define rssringoccs_Sqrt_Pi_By_Two          1.25331413731550025
#define rssringoccs_Sqrt_One_By_Pi          0.56418958354775628
#define rssringoccs_Sqrt_Two_By_Pi          0.79788456080286535
#define rssringoccs_Two_By_Sqrt_Pi          1.12837916709551257
#define rssringoccs_Pi_By_Two               1.57079632679489661
#define rssringoccs_Three_Pi_By_Four        2.35619449019234492
#define rssringoccs_Pi_By_Four              0.78539816339744830
#define rssringoccs_One_Pi                  3.14159265358979323
#define rssringoccs_Two_Pi                  6.28318530717958647
#define rssringoccs_Sqrt_Two                1.41421356237309504
#define rssringoccs_Rcpr_Euler_E            0.36787944117144232
#define rssringoccs_Euler_E                 2.71828182845904523
#define rssringoccs_Natural_Log_of_2        0.69314718055994531
#define rssringoccs_Natural_Log_Of_10       2.30258509299404568

/*  Long double precision constants.                                          */
#define rssringoccs_Sqrt_One_By_Two_Pi_L    0.3989422804014326779399461L
#define rssringoccs_Sqrt_Pi_By_Eight_L      0.6266570686577501256039413L
#define rssringoccs_Sqrt_Pi_By_Two_L        1.2533141373155002512078830L
#define rssringoccs_Sqrt_One_By_Pi_L        0.5641895835477562869480795L
#define rssringoccs_Sqrt_Two_By_Pi_L        0.7978845608028653558798921L
#define rssringoccs_Two_By_Sqrt_Pi_L        1.1283791670955125738961590L
#define rssringoccs_Pi_By_Two_L             1.5707963267948966192313220L
#define rssringoccs_Three_Pi_By_Four_L      2.3561944901923449288469830L
#define rssringoccs_Pi_By_Four_L            0.7853981633974483096156608L
#define rssringoccs_One_Pi_L                3.1415926535897932384626430L
#define rssringoccs_Two_Pi_L                6.2831853071795864769252870L
#define rssringoccs_Sqrt_Two_L              1.4142135623730950488016890L
#define rssringoccs_Rcpr_Euler_E_L          0.3678794411714423215955238L
#define rssringoccs_Euler_E_L               2.7182818284590452353602875L
#define rssringoccs_Natural_Log_of_2_L      0.6931471805599453094172321L
#define rssringoccs_Natural_Log_Of_10_L     2.3025850929940456840179910L

/*  Degrees to radians conversion, and vice-versa.                            */
#define rssringoccs_Deg_To_Rad_F 1.745329251e-2F
#define rssringoccs_Deg_To_Rad   1.74532925199432957e-2
#define rssringoccs_Deg_To_Rad_L 1.74532925199432957692369076e-2L

#define rssringoccs_Rad_To_Deg_F 57.295779513F
#define rssringoccs_Rad_To_Deg   57.29577951308232087
#define rssringoccs_Rad_To_Deg_L 57.29577951308232087679815481L

/*  The speed of light in km/s.                                               */
#define rssringoccs_Speed_Of_Light_KMS_F  299792.4580F
#define rssringoccs_Speed_Of_Light_KMS    299792.4580
#define rssringoccs_Speed_Of_Light_KMS_L  299792.4580L

/*  Macros for the largest values of float, double, and long double,          *
 *  respectively, that will not return INFINITY when exp(x) is computed.      */
#define rssringoccs_Max_Float_Base_E    (FLT_MAX_10_EXP *                      \
                                         rssringoccs_Natural_Log_Of_10_F)

#define rssringoccs_Max_Double_Base_E   (DBL_MAX_10_EXP *                      \
                                         rssringoccs_Natural_Log_Of_10)

#define rssringoccs_Max_LDouble_Base_E  (LDBL_MAX_10_EXP *                     \
                                         rssringoccs_Natural_Log_Of_10_L)

/*  Macros for the smallest values of float, double, and long double,         *
 *  respectively, that will not return 0 when exp(x) is computed.             */
#define rssringoccs_Min_Float_Base_E    (FLT_MIN_10_EXP *                      \
                                         rssringoccs_Natural_Log_Of_10_F)

#define rssringoccs_Min_Double_Base_E   (DBL_MIN_10_EXP *                      \
                                         rssringoccs_Natural_Log_Of_10)

#define rssringoccs_Min_LDouble_Base_E  (LDBL_MIN_10_EXP *                     \
                                         rssringoccs_Natural_Log_Of_10_L)

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Double_Abs                                                *
 *  Purpose:                                                                  *
 *      Compute the absolute value of a real number (fabs alias).             *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          A real number.                                                    *
 *  Output:                                                                   *
 *      double abs_x:                                                         *
 *          The absolute value of x, |x|.                                     *
 ******************************************************************************/
extern float rssringoccs_Float_Abs(float x);
extern double rssringoccs_Double_Abs(double x);
extern long double rssringoccs_LDouble_Abs(long double x);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Double_Arctan                                             *
 *  Purpose:                                                                  *
 *      Compute the arctan (inverse tangent) of a real number.                *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          A real number.                                                    *
 *  Output:                                                                   *
 *      double atan_x:                                                        *
 *          The inverse tangent of x, tan^-1(x).                              *
 ******************************************************************************/
extern float rssringoccs_Float_Arctan(float x);
extern double rssringoccs_Double_Arctan(double x);
extern long double rssringoccs_LDouble_Arctan(long double x);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Double_Arctan2                                            *
 *  Purpose:                                                                  *
 *      Compute the 2-dimensional arctan (inverse tangent) of a point in the  *
 *      Cartesian plane. This is the angle the point makes with the positive  *
 *      x-axis.                                                               *
 *  Arguments:                                                                *
 *      double y:                                                             *
 *          A real number.                                                    *
 *      double x:                                                             *
 *          Another real number.                                              *
 *  Output:                                                                   *
 *      double atan:                                                          *
 *          The angle the point (x,y) makes with (1,0) in the plane.          *
 *  NOTES:                                                                    *
 *      By convention dating back to (at least) the 1970s, Arctan2 takes the  *
 *      input as (y,x), not (x,y). i.e. the first argument is the y           *
 *      component and the second argument is the x component. This is contrary*
 *      to most 2 dimensional functions that want their inputs as (x,y).      *
 *      This is probably because we are trying to compute tan^-1(y/x) but     *
 *      need to be careful about the signs of y and x, so we write            *
 *      arctan(y,x).                                                          *
 *                                                                            *
 *      This returns a number between -pi and pi, so there is a "branch cut"  *
 *      along the negative x axis. Because of this, use of this function      *
 *      in complex routines results in actual branch cuts.                    *
 ******************************************************************************/
extern float rssringoccs_Float_Arctan2(float y, float x);
extern double rssringoccs_Double_Arctan2(double y, double x);
extern long double rssringoccs_LDouble_Arctan2(long double y, long double x);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Double_Copysign                                           *
 *  Purpose:                                                                  *
 *      Given two numbers x and y, returns a value that has the magnitude of  *
 *      x and the sign of y.                                                  *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          A real number.                                                    *
 *      double y:                                                             *
 *          Another real number.                                              *
 *  Output:                                                                   *
 *      double z:                                                             *
 *          The value sgn(y) * |x|.                                           *
 ******************************************************************************/
extern float rssringoccs_Float_Copysign(float x, float y);
extern double rssringoccs_Double_Copysign(double x, double y);
extern long double rssringoccs_LDouble_Copysign(long double x, long double y);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Double_Cos                                                *
 *  Purpose:                                                                  *
 *      Computes the cosine of a real number.                                 *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          A real number.                                                    *
 *  Output:                                                                   *
 *      double cos_x:                                                         *
 *          The cosine of x, cos(x).                                          *
 ******************************************************************************/
extern float rssringoccs_Float_Cos(float x);
extern double rssringoccs_Double_Cos(double x);
extern long double rssringoccs_LDouble_Cos(long double x);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Double_Cosh                                               *
 *  Purpose:                                                                  *
 *      Computes the hyperbolic cosine of a real number.                      *
 *  Arguments:                                                                *
 *      double x:                                                             *
 *          A real number.                                                    *
 *  Output:                                                                   *
 *      double cos_x:                                                         *
 *          The hyperbolic cosine of x, cosh(x).                              *
 ******************************************************************************/
extern float rssringoccs_Float_Cosh(float x);
extern double rssringoccs_Double_Cosh(double x);
extern long double rssringoccs_LDouble_Cosh(long double x);

/*  Aliases for the sine trig function found in math.h.                       */
extern float rssringoccs_Float_Sin(float x);
extern double rssringoccs_Double_Sin(double x);
extern long double rssringoccs_LDouble_Sin(long double x);

/*  Aliases for the cosine tan function found in math.h.                      */
extern float rssringoccs_Float_Tan(float x);
extern double rssringoccs_Double_Tan(double x);
extern long double rssringoccs_LDouble_Tan(long double x);

/*  Aliases for the square root function found in math.h.                     */
extern float rssringoccs_Float_Sqrt(float x);
extern double rssringoccs_Double_Sqrt(double x);
extern long double rssringoccs_LDouble_Sqrt(long double x);

/*  Aliases for the exponential function found in math.h.                     */
extern float rssringoccs_Float_Exp(float x);
extern double rssringoccs_Double_Exp(double x);
extern long double rssringoccs_LDouble_Exp(long double x);

/*  Aliases for the exponential function found in math.h.                     */
extern float rssringoccs_Float_Log(float x);
extern double rssringoccs_Double_Log(double x);
extern long double rssringoccs_LDouble_Log(long double x);

/*  Set INFINITY to the HUGE_VAL macro that is specified in math.h. Most      *
 *  implementations already have an INFINITY macro, but it is not required.   */
#define rssringoccs_Infinity (HUGE_VAL)
#define rssringoccs_Infinity_F ((float)(rssringoccs_Infinity))
#define rssringoccs_Infinity_L ((long double)(rssringoccs_Infinity))

/*  We'll use the CPYTHON method of defining NAN, the source code of which is *
 *   contained in python/cpython/Include/pymath.h.                            */
#define rssringoccs_NaN (rssringoccs_Infinity * 0.0)
#define rssringoccs_NaN_F ((float)(rssringoccs_NaN))
#define rssringoccs_NaN_L ((long double)(rssringoccs_NaN))

#define rssringoccs_Is_Inf(x) ((x) == ((x)+1))
#define rssringoccs_Is_NaN(x) ((x) != (x))

/*  The following functions are not required in C89/C90 implementations of    *
 *  math.h. The algorithms for their computations are very straight-forward,  *
 *  reducing the definitions to computations using sin, cos, exp, etc. We     *
 *  provide them here for portability.                                        */

extern float rssringoccs_Float_Sinc(float x);
extern double rssringoccs_Double_Sinc(double x);
extern long double rssringoccs_LDouble_Sinc(long double x);

extern float rssringoccs_Float_Sinh(float x);
extern double rssringoccs_Double_Sinh(double x);
extern long double rssringoccs_LDouble_Sinh(long double x);

extern float rssringoccs_Float_Tanh(float x);
extern double rssringoccs_Double_Tanh(double x);
extern long double rssringoccs_LDouble_Tanh(long double x);

extern float rssringoccs_Float_Erf(float x);
extern double rssringoccs_Double_Erf(double x);
extern long double rssringoccs_LDouble_Erf(long double x);

extern float rssringoccs_Float_Erfc(float x);
extern double rssringoccs_Double_Erfc(double x);
extern long double rssringoccs_LDouble_Erfc(long double x);

extern float rssringoccs_Float_Erfcx(float x);
extern double rssringoccs_Double_Erfcx(double x);
extern long double rssringoccs_LDouble_Erfcx(long double x);

extern float rssringoccs_Float_Faddeeva_Im(float x);
extern double rssringoccs_Double_Faddeeva_Im(double x);
extern long double rssringoccs_LDouble_Faddeeva_Im(long double x);

extern unsigned long rssringoccs_Factorial(unsigned int n);

extern unsigned long
rssringoccs_Falling_Factorial(unsigned int x, unsigned int N);

extern float
rssringoccs_Real_Poly_Float_Coeffs(float *coeffs, unsigned int degree, float x);

extern double
rssringoccs_Real_Poly_Double_Coeffs(double *coeffs,
                                    unsigned int degree,
                                    double x);

extern long double
rssringoccs_Real_Poly_LDouble_Coeffs(long double *coeffs,
                                     unsigned int degree,
                                     long double x);

extern float
rssringoccs_Real_Poly_Deriv_Float_Coeffs(float *coeffs, unsigned int degree,
                                         unsigned int deriv, float x);

extern double
rssringoccs_Real_Poly_Deriv_Double_Coeffs(double *coeffs, unsigned int degree,
                                          unsigned int deriv, double x);

extern long double
rssringoccs_Real_Poly_Deriv_LDouble_Coeffs(long double *coeffs,
                                           unsigned int degree,
                                           unsigned int deriv,
                                           long double x);

extern void
rssringoccs_Reverse_Float_Array(float *arr, unsigned long arrsize);

extern void
rssringoccs_Reverse_Double_Array(double *arr, unsigned long arrsize);

extern void
rssringoccs_Reverse_LDouble_Array(long double *arr, unsigned long arrsize);

#endif
/*  End of include guard.                                                     */
