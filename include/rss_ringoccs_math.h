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
 *      This header file allows for compatibility of rss_ringoccs for users   *
 *      of the C89/C90 math.h header file, and the C99/C11 math.h. The later  *
 *      versions of math.h provide the following:                             *
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
 *      This file also provides NAN and INFINITY macros if they are not set.  *
 *      NOTE:                                                                 *
 *          INFINITY is set as the standard macro HUGE_VAL defined in math.h  *
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
 *          Similarly, NAN is defined as HUGE_VAL * 0, which should be        *
 *          infinity times zero, which is Not-A-Number. Python does this as   *
 *          well, issuing the following warning:                              *
 *              Py_NAN                                                        *
 *              A value that evaluates to a NaN. On IEEE 754 platforms INF*0  *
 *              or INF/INF works. Define Py_NO_NAN in pyconfig.h if your      *
 *              platform doesn't support NaNs.                                *
 *          If necessary, redefine NAN here to whatever your platform allows. *
 *      Lastly, this file provides a bunch of constants that are commonly     *
 *      used, as well as various math functions that are not included in      *
 *      either C89 or C99 math.h.                                             *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) math.h:                                                               *
 *      Standard library for mathematical functions like sin, cos, atan.      *
 *  2.) float.h:                                                              *
 *      Standard library which contains macros for the smallest and largest   *
 *      values allowed by your system.                                        *
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

/*  Include guard for this file to prevent including this twice.              */
#ifndef __RSS_RINGOCCS_MATH_H__
#define __RSS_RINGOCCS_MATH_H__

/*  The __RSS_RINGOCCS_USING_C99_MATH_H__ macro is found here.                */
#include <rss_ringoccs/include/rss_ringoccs_config.h>

/*  It is not assumed you have C99 math.h, but rather the C89 math.h, which   *
 *  is also referred to as the ANSI C version. If you have C99 available,     *
 *  change __HAS_C99_MATH_H__ in rss_ringoccs_config.h to 1 and then          *
 *  rebuild rss_ringoccs. Some implementations of the C99 standard don't      *
 *  define float/long double versions of math functions and simply do         *
 *  something like sinf(x) = (float)sin((double)x), so there is no difference.*/

/*  Check which version of C you are compiling with and set the macro         *
 *  __HAS_C99_MATH_H__ accordingly.                                           */
#if __RSS_RINGOCCS_USING_C99_MATH_H__ == 1

#if __STDC_VERSION__ >= 199901L
#define __HAS_C99_MATH_H__ 1

#else
/*  Else statement for #if __STDC_VERSION__ >= 199901L.                       */

/*  You requested C99 math.h but don't have it. Abort with error.             */
#error Requested C99 math.h but your compiler does not support it.

#endif
/*  End of #if __STDC_VERSION__ >= 199901L                                    */

#else
/*  Else statement for #if __RSS_RINGOCCS_USING_C99_MATH_H__ == 1.            */

#define __HAS_C99_MATH_H__ 0

#endif
/*  End of #if __RSS_RINGOCCS_USING_C99_MATH_H__ == 1.                        */

/*  Include the standard library header math.h. We're only going to alias     *
 *  functions we ever use in rss_ringoccs, sin, cos, fabs, exp, atan2.        */
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
#define rssringoccs_Py_By_Four_F            0.785398163F
#define rssringoccs_One_Pi_F                3.141592654F
#define rssringoccs_Two_Pi_F                6.283185307F
#define rssringoccs_Sqrt_Two_F              1.414213562F
#define rssringoccs_Rcpr_Euler_E_F          0.367879441F
#define rssringoccs_Euler_E_F               2.718281828F
#define rssringoccs_Natural_Log_Of_10_F     2.302585093F

/*  Double precision constants.                                               */
#define rssringoccs_Sqrt_One_By_Two_Pi      0.39894228040143267
#define rssringoccs_Sqrt_Pi_By_Eight        0.62665706865775012
#define rssringoccs_Sqrt_Pi_By_Two          1.25331413731550025
#define rssringoccs_Sqrt_One_By_Pi          0.56418958354775628
#define rssringoccs_Sqrt_Two_By_Pi          0.79788456080286535
#define rssringoccs_Two_By_Sqrt_Pi          1.12837916709551257
#define rssringoccs_Pi_By_Two               1.57079632679489661
#define rssringoccs_Py_By_Four              0.78539816339744830
#define rssringoccs_One_Pi                  3.14159265358979323
#define rssringoccs_Two_Pi                  6.28318530717958647
#define rssringoccs_Sqrt_Two                1.41421356237309504
#define rssringoccs_Rcpr_Euler_E            0.36787944117144232
#define rssringoccs_Euler_E                 2.71828182845904523
#define rssringoccs_Natural_Log_Of_10       2.30258509299404568

/*  Long double precision constants.                                          */
#define rssringoccs_Sqrt_One_By_Two_Pi_L    0.3989422804014326779399461L
#define rssringoccs_Sqrt_Pi_By_Eight_L      0.6266570686577501256039413L
#define rssringoccs_Sqrt_Pi_By_Two_L        1.2533141373155002512078830L
#define rssringoccs_Sqrt_One_By_Pi_L        0.5641895835477562869480795L
#define rssringoccs_Sqrt_Two_By_Pi_L        0.7978845608028653558798921L
#define rssringoccs_Two_By_Sqrt_Pi_L        1.1283791670955125738961590L
#define rssringoccs_Pi_By_Two_L             1.5707963267948966192313220L
#define rssringoccs_Py_By_Four_L            0.7853981633974483096156608L
#define rssringoccs_One_Pi_L                3.1415926535897932384626430L
#define rssringoccs_Two_Pi_L                6.2831853071795864769252870L
#define rssringoccs_Sqrt_Two_L              1.4142135623730950488016890L
#define rssringoccs_Rcpr_Euler_E_L          0.3678794411714423215955238L
#define rssringoccs_Euler_E_L               2.7182818284590452353602875L
#define rssringoccs_Natural_Log_Of_10_L     2.3025850929940456840179910L

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

/*  Aliases for the sine trig function found in math.h.                       */
extern float rssringoccs_Float_Sin(float x);
extern double rssringoccs_Double_Sin(double x);
extern long double rssringoccs_LDouble_Sin(long double x);

/*  Aliases for the cosine trig function found in math.h.                     */
extern float rssringoccs_Float_Cos(float x);
extern double rssringoccs_Double_Cos(double x);
extern long double rssringoccs_LDouble_Cos(long double x);

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

/*  Aliases for the absolute value function found in math.h.                  */
extern float rssringoccs_Float_Abs(float x);
extern double rssringoccs_Double_Abs(double x);
extern long double rssringoccs_LDouble_Abs(long double x);

/*  Aliases for the atan function found in math.h.                            */
extern float rssringoccs_Float_Arctan(float x);
extern double rssringoccs_Double_Arctan(double x);
extern long double rssringoccs_LDouble_Arctan(long double x);

/*  Aliases for the atan2 function found in math.h.                           */
extern float rssringoccs_Float_Arctan2(float y, float x);
extern double rssringoccs_Double_Arctan2(double y, double x);
extern long double rssringoccs_LDouble_Arctan2(long double y, long double x);

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

extern float rssringoccs_Float_Sinh(float x);
extern double rssringoccs_Double_Sinh(double x);
extern long double rssringoccs_LDouble_Sinh(long double x);

extern float rssringoccs_Float_Cosh(float x);
extern double rssringoccs_Double_Cosh(double x);
extern long double rssringoccs_LDouble_Cosh(long double x);

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

extern double rssringoccs_Double_Copysign(double x, double y);

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

#endif
/*  End of include guard.                                                     */
