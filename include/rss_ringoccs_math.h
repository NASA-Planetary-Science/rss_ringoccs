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
 *      (it's defined in this file before the #include <math.h> statement).   *
 *      The functions have the following aliases if __HAS_C99_MATH_H__ is set:*
 *          Func_Float      = funcf                                           *
 *          Func_Double     = func                                            *
 *          Func_LongDouble = funcl                                           *
 *      If not, we do type conversion:                                        *
 *          Func_Float(x)      = func((double)(x))                            *
 *          Func_Double(x)     = func((x))                                    *
 *          Func_LongDouble(x) = func((double)(x))                            *
 *      Note the extra parentheses surrounding the "x" in the func. This is   *
 *      because the preprocessor takes whatever you give it and write this    *
 *      verbatim into the expression provided. If we did:                     *
 *          Func_Float(x) = func((double)x)                                   *
 *      Then doing something like Func_Float(x+y) will return                 *
 *      func((double)x+y), so only x gets converted to double, which defeats  *
 *      the purpose. With the extra parentheses we get                        *
 *      Func_Float(x+y) = func((double)(x+y)), which is what we want.         *
 *                                                                            *
 *      This also provides NAN and INFINITY macros if they are not set.       *
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
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  If your implementation supports the C89/C90 math.h header file, and is    *
 *  lacking the extra function declarations, you will need to comment out the *
 *  #define __HAS_C99_MATH_H__ statement or you will have compiler errors.    *
 *  If HUGE_VAL and HUGE_VAL*0 don't correspond to infinity and nan,          *
 *  respectively, on your system you will need to alter those macros. This    *
 *  has worked on macOS and linux systems when tested, however the GNU glibc  *
 *  implementation provides those macros already, so no extra effort needed.  *
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
#ifndef _RSS_RINGOCCS_MATH_H_
#define _RSS_RINGOCCS_MATH_H_

/*  It is not assumed you have C99 math.h, but rather the C89 math.h, which   *
 *  is also (incorrectly) referred to as the ANSI C version. If you have C99  *
 *  available, change the __HAS_C99_MATH_H__ to 1 and rebuild rss_ringoccs.   */

/*  Check which version of C you are compiling with and set the macro         *
 *  __HAS_C99_MATH_H__ accordingly.                                           */
#if __STDC_VERSION__ >= 199901L
#define __HAS_C99_MATH_H__ 1
#else
#define __HAS_C99_MATH_H__ 0
#endif
/*  End of #if __STDC_VERSION__ >= 199901L                                    */

/*  Include the standard library header math.h. We're only going to alias     *
 *  functions we ever use in rss_ringoccs, sin, cos, fabs, exp, atan2.        */
#include <math.h>

/*  And this header file contains macros for the smallest and largest allowed *
 *  values for your system.                                                   */
#include <float.h>

/*  This number is needed for computing the smallest value which will yield   *
 *  INFINITY or HUGE_VAL with the exponential function.                       */
#define NATURAL_LOG_OF_10 2.302585092994045684017991

/* Define Miscellaneous Constants.                                            */
#define SQRT_ONE_BY_TWO_PI  0.3989422804014326779399461
#define SQRT_PI_BY_8        0.6266570686577501256039413
#define SQRT_PI_BY_2        1.2533141373155002512078830
#define SQRT_ONE_BY_PI      0.5641895835477562869480794
#define SQRT_TWO_BY_PI      0.7978845608028653558798921
#define PI_BY_TWO           1.570796326794896619231321691639751440
#define PI_BY_FOUR          0.785398163397448309615660845819875721
#define ONE_PI              3.141592653589793238462643383279502880
#define TWO_PI              6.283185307179586476925286766559005768
#define SQRT_2              1.414213562373095048801688724209698080
#define RCPR_EULER_E        0.367879441171442321595523770161460867
#define ONE_BY_SQRT_PI      0.564189583547756286948079451560772585
#define TWO_BY_SQRT_PI      1.128379167095512573896158903121545171

/*  Define the speed of light in kilometers per second.                       */
#define SPEED_OF_LIGHT_KMS  299792.4580

/*  Macros for the largest values of float, double, and long double,          *
 *  respectively, that will not return INFINITY when exp(x) is computed.      */
#define MAX_FLOAT_BASE_E    FLT_MAX_10_EXP * NATURAL_LOG_OF_10
#define MAX_DOUBLE_BASE_E   DBL_MAX_10_EXP * NATURAL_LOG_OF_10
#define MAX_LDOUBLE_BASE_E  LDBL_MAX_10_EXP * NATURAL_LOG_OF_10

/*  Aliases for the sine trig function found in math.h.                       */
extern float rssringoccs_Float_Sin(float x);
extern double rssringoccs_Double_Sin(double x);
extern long double rssringoccs_LongDouble_Sin(long double x);

/*  Aliases for the cosine trig function found in math.h.                     */
extern float rssringoccs_Float_Cos(float x);
extern double rssringoccs_Double_Cos(double x);
extern long double rssringoccs_LongDouble_Cos(long double x);

/*  Aliases for the cosine tan function found in math.h.                      */
#define rssringoccs_Tan_Double(x)       tan(x)
#if __HAS_C99_MATH_H__ != 0
#define rssringoccs_Tan_Float(x)        tanf(x)
#define rssringoccs_Tan_LongDouble(x)   tanl(x)
#else
#define rssringoccs_Tan_Float(x)        tan((double)(x))
#define rssringoccs_Tan_LongDouble(x)   tan((double)(x))
#endif

/*  Aliases for the square root function found in math.h.                     */
extern float rssringoccs_Float_Sqrt(float x);
extern double rssringoccs_Double_Sqrt(double x);
extern long double rssringoccs_LongDouble_Sqrt(long double x);

/*  Aliases for the exponential function found in math.h.                     */
extern float rssringoccs_Float_Exp(float x);
extern double rssringoccs_Double_Exp(double x);
extern long double rssringoccs_LongDouble_Exp(long double x);

/*  Aliases for the exponential function found in math.h.                     */
extern float rssringoccs_Float_Log(float x);
extern double rssringoccs_Double_Log(double x);
extern long double rssringoccs_LongDouble_Log(long double x);

/*  Aliases for the absolute value function found in math.h.                  */
extern float rssringoccs_Float_Abs(float x);
extern double rssringoccs_Double_Abs(double x);
extern long double rssringoccs_LongDouble_Abs(long double x);

/*  Aliases for the atan function found in math.h.                            */
extern float rssringoccs_Float_Arctan(float x);
extern double rssringoccs_Double_Arctan(double x);
extern long double rssringoccs_LongDouble_Arctan(long double x);

/*  Aliases for the atan2 function found in math.h.                           */
extern float rssringoccs_Float_Arctan2(float y, float x);
extern double rssringoccs_Double_Arctan2(double y, double x);
extern long double rssringoccs_LongDouble_Arctan2(long double y, long double x);

/*  If INFINITY is not defined, set it to the HUGE_VAL macro that is          *
 *  specified in math.h. Most compilers already have an INFINITY macro, but   *
 *  it is not required in the C89 standard.                                   */
#ifndef INFINITY
#define rssringoccs_Infinity HUGE_VAL
#else
#define rssringoccs_Infinity INFINITY
#endif

/*  If NAN is not defined, we'll use the CPYTHON method of defining NAN, the  *
 *  source code of which is contained in python/cpython/Include/pymath.h.     */
#ifndef NAN
#define rssringoccs_NaN (INFINITY * 0.0)
#else
#define rssringoccs_NaN NAN
#endif

/*  The following functions are not required in C89/C90 implementations of    *
 *  math.h. The algorithms for their computations are very straight-forward,  *
 *  reducing the definitions to computations using sin, cos, exp, etc. We     *
 *  provide them here for portability.                                        */

extern float rssringoccs_Float_Sinh(float x);
extern double rssringoccs_Double_Sinh(double x);
extern long double rssringoccs_LongDouble_Sinh(long double x);

extern float rssringoccs_Float_Cosh(float x);
extern double rssringoccs_Double_Cosh(double x);
extern long double rssringoccs_LongDouble_Cosh(long double x);

extern float rssringoccs_Float_Tanh(float x);
extern double rssringoccs_Double_Tanh(double x);
extern long double rssringoccs_LongDouble_Tanh(long double x);

extern float rssringoccs_Float_Erf(float x);
extern double rssringoccs_Double_Erf(double x);
extern long double rssringoccs_LongDouble_Erf(long double x);

extern float rssringoccs_Float_Erfc(float x);
extern double rssringoccs_Double_Erfc(double x);
extern long double rssringoccs_LongDouble_Erfc(long double x);

extern float rssringoccs_Float_Erfcx(float x);
extern double rssringoccs_Double_Erfcx(double x);
extern long double rssringoccs_LongDouble_Erfcx(long double x);

extern float rssringoccs_Float_Faddeeva_Im(float x);
extern double rssringoccs_Double_Faddeeva_Im(double x);
extern long double rssringoccs_LongDouble_Faddeeva_Im(long double x);

extern long rssringoccs_Factorial(int n);
extern long rssringoccs_Falling_Factorial(int x, int N);

extern float
rssringoccs_Float_Poly_Float_Coeffs(float *coeffs,
                                    unsigned int degree, float x);

extern double
rssringoccs_Double_Poly_Double_Coeffs(double *coeffs,
                                      unsigned int degree, double x);

extern long double
rssringoccs_LongDouble_Poly_LongDouble_Coeffs(long double *coeffs,
                                              unsigned int degree,
                                              long double x);

extern float
rssringoccs_Float_Poly_Deriv_Float_Coeffs(float *coeffs,
                                          unsigned int degree,
                                          unsigned int deriv,
                                          float x);

extern double
rssringoccs_Double_Poly_Deriv_Double_Coeffs(double *coeffs,
                                            unsigned int degree,
                                            unsigned int deriv,
                                            double x);

extern long double
rssringoccs_LongDouble_Poly_Deriv_LongDouble_Coeffs(long double *coeffs,
                                                    unsigned int degree,
                                                    unsigned int deriv,
                                                    long double x);


#endif
/*  End of include guard.                                                     */
