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
 *                       rss_ringoccs_complex_erfc                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complementary complex error function.    *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CFloat_Erfc:                                              *
 *      rssringoccs_CDouble_Erfc:                                             *
 *      rssringoccs_CLDouble_Erfc:                                            *
 *  Purpose:                                                                  *
 *      Computes the complementary error function of a complex value z.       *
 *                                                                            *
 *                            infinity                                        *
 *                               -                                            *
 *                      2       | |                                           *
 *          Erf(z) = -------    |   exp(-t^2) dt = 1 - Erf(z)                 *
 *                   sqrt(pi) | |                                             *
 *                             -                                              *
 *                             z                                              *
 *                                                                            *
 *      Note, the integral formula only applied to real z. For complex we use *
 *      the expression 1 - Erf(z) as the definition.                          *
 *  Arguments:                                                                *
 *      z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      erfc_z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):    *
 *          The complementary error function, Erfc(z).                        *
 *  Method:                                                                   *
 *      Analyze the input to see which region of the complex plane it lies in *
 *      and use the fadeeva and Erfcx functions to compute.                   *
 *  NOTES:                                                                    *
 *      No actual float or long double algorithms have been implemented by    *
 *      rss_ringoccs. The complementary error functions simply convert float  *
 *      and long double inputs to doubles.                                    *
 *                                                                            *
 *      This is an alteration of the MIT Faddeeva package. It has been        *
 *      altered to be C89/C90 compliant and uses the rest of the rss_ringoccs *
 *      complex routines to perform the computation. The original authors are:*
 *          Steven G. Johnson, Massachusetts Institute of Technology          *
 *          Joachim Wuttke, Forschungszentrum Jülich, 2013                    *
 *      Original license and copyright notice found below. Most of the        *
 *      original code has been rewritten so that complex.h/C99 features are   *
 *      NOT required. This file compiles with gcc using the -std=c89 -ansi    *
 *      -Wall -Wextra and -Wpedantic options.                                 *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 *  2.) rss_ringoccs_math.h:                                                  *
 *          Header file containing lots of real-valued math functions.        *
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
 *  Date:       December 23, 2020                                             *
 ******************************************************************************/

/******************************************************************************
 *                             LIBCERF LICENSE                                *
 ******************************************************************************/

/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File cerf.h:
 *   Declare exported functions.
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   ../CHANGELOG
 *
 * Man pages:
 *   w_of_z(3), dawson(3), voigt(3), cerf(3), erfcx(3), erfi(3)
 */


/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Double precision complementary error function.                            */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Erfc(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double z_x, z_y, w_x, w_y;
    double mRe_z2, mIm_z2;
    rssringoccs_ComplexDouble w, temp;

    /*  Extract the real and imaginary parts from the input z.                */
    z_x = rssringoccs_CDouble_Real_Part(z);
    z_y = rssringoccs_CDouble_Imag_Part(z);

    /*  For purely imaginary inputs we can use the Faddeeva function. This    *
     *  is a scaled error function defined by w(z) = exp(-z^2) Erfc(-iz).     */
    if (z_x == 0.0)
    {
        /*  Set the real part to 1. This is because Erfc(z) = 1 - Erf(z) and  *
         *  Erf(iy) has zero real part. Hence the real part of Erfc(iy) is 1. */
        w_x = 1.0;

        /*  For large values of the imaginary part, Erfc tends to infinity.   *
         *  Since the imaginary part of the Faddeeva function tends to zero   *
         *  we'd end up with infinity*0 resulting in a NaN. We handle this    *
         *  case seperately.                                                  */

        /*  If the square of the y component is greater than the              *
         *  rssringoccs_Max_Double_Base_E macro, exp(y^2) will return         *
         *  infinity. Check this.                                             */
        if (z_y*z_y > rssringoccs_Max_Double_Base_E)
        {
            /*  Erf(0+iy) tends to infinity for large y, so                   *
             *  Erfc(0+iy) = 1 - Erf(0+iy) tends to 1 - i infinity.           */
            if (z_y>0)
                w_y = -rssringoccs_Infinity;

            /*  The error function is odd valued, so Erf(0 - iy) tends to     *
             *  minus infinity as y gets large and positive. Hence Erfc tends *
             *  to 1 + i infinity.                                            */
            else
                w_y = rssringoccs_Infinity;
        }

        /*  If exp(y^2) won't return infinity, use the definition of the      *
         *  Faddeeva function to compute the imaginary part of Erfc(z).       */
        else
            w_y = -exp(z_y*z_y) * rssringoccs_Double_Faddeeva_Im(z_y);

        /*  Create the output complex number using rssringoccs_CDouble_Rect.  */
        w = rssringoccs_CDouble_Rect(w_x, w_y);
    }
    /*  End of if (z_x == 0.0).                                               */

    /*  If the input is purely real, we can use the real valued complementary *
     *  error function, which can be computed from the Erfcx function. Erfcx  *
     *  is a scaled error function defined by erfcx(x) = exp(x^2)*erfc(x).    */
    else if (z_y == 0.0)
    {
        /*  The imaginary part of erfc(x) is also zero since erf(x) has no    *
         *  imaginary part and erfc(x) = 1-erf(x). Set w_y to minus z_y in    *
         *  case a "signed zero" is needed in a computation.                  */
        w_y = -z_y;

        /*  Check for underflow. The macro rssringoccs_Min_Double_Base_E is   *
         *  the smallest value which will not yield zero when exp(x) is       *
         *  computed. We need to check if -x^2 will result in and underflow.  */
        if (-z_x*z_x < rssringoccs_Min_Double_Base_E)
        {
            /*  If x is positive and -x^2 < rssringoccs_Min_Double_Base_E,    *
             *  then erf(x) is 1 to many decimals, and erfc(x) = 0.           */
            if (z_x >= 0.0)
                w_x = 0.0;

            /*  For negative x erf(x) is roughly -1, so erfc(x) = 1-(-1) = 2. */
            else
                w_x = 2.0;
        }

        /*  If x^2 is small enough that we won't underflow, use the Erfcx     *
         *  function and apply the formula described above.                   */
        else
        {
            /*  For positive x we can use erfc(x) = exp(-x^2) erfcx(x).       */
            if (z_x >= 0.0)
                w_x = exp(-z_x*z_x) * rssringoccs_Double_Erfcx(z_x);

            /*  And for negative we use the reflection formula for erfc.      */
            else
                w_x = 2.0 - exp(-z_x*z_x) * rssringoccs_Double_Erfcx(z_x);
        }

        /*  Create the output complex number using rssringoccs_CDouble_Rect.  */
        w = rssringoccs_CDouble_Rect(w_x, w_y);
    }
    /*  End of else statement else if (z_y == 0.0).                           */

    /*  For the more general z = x + iy where x and y are non-zero we'll use  *
     *  the complex Faddeeva function w(z). This is defined by                *
     *  w(z) = exp(-z^2) erfc(-iz).                                           */
    else
    {
        /*  Compute the real part of -z^2 = -(x+iy)^2 = y^2 - x^2 - i2xy. We  *
         *  can save a multiplication by writing y^2-x^2 = (y-x)(x+y). This   *
         *  also has the added caution of potentially avoiding underflow.     */
        mRe_z2 = (z_y - z_x) * (z_x + z_y);

        /*  Similarly, compute the imaginary part of -z^2.                    */
        mIm_z2 = -2.0*z_x*z_y;

        /*  In the underflow case where the real part has a value less then   *
         *  the macro rssringoccs_Min_Double_Base_E we'll simply return       *
         *  2.0. We have erfc(z) = exp(-z^2) w(iz). Expanding this yields:    *
         *                                                                    *
         *      erfc(z) = exp(-mRe_z2 + i mIm_z2) w(iz)                       *
         *              = exp(-mRe_z2)(cos(mIm_z2) + i sin(mIm_z2)) w(iz)     *
         *                                                                    *
         *  The exp term will dominate the imaginary part and the real part   *
         *  is determined by the sign of x.                                   */
        if (mRe_z2 < rssringoccs_Min_Double_Base_E)
        {
            w_y = 0.0;
            if (z_x >= 0.0)
                w_x = 0.0;
            else
                w_x = 2.0;
            w = rssringoccs_CDouble_Rect(w_x, w_y);
        }

        /*  If we won't underflow, use the complex Faddeeva function.         */
        else
        {
            if (z_x >= 0.0)
            {
                temp = rssringoccs_CDouble_Rect(mRe_z2, mIm_z2);
                temp = rssringoccs_CDouble_Exp(temp);
                w = rssringoccs_CDouble_Rect(-z_y, z_x);
                w = rssringoccs_CDouble_Faddeeva(w);
                w = rssringoccs_CDouble_Multiply(w, temp);
            }

            /*  For negative x, use the reflection formula.                   */
            else
            {
                temp = rssringoccs_CDouble_Rect(mRe_z2, mIm_z2);
                temp = rssringoccs_CDouble_Exp(temp);
                w = rssringoccs_CDouble_Rect(z_y, -z_x);
                w = rssringoccs_CDouble_Faddeeva(w);
                w = rssringoccs_CDouble_Multiply(w, temp);
                w = rssringoccs_CDouble_Subtract_Real(2.0, w);
            }
        }
        /*  End of if (mRe_z2 < rssringoccs_Min_Double_Base_E).               */
    }
    /*  End of else statement for if (z_x == 0.0).                            */

    return w;
}
/*  End of rssringoccs_CDouble_Erfc.                                          */

/*  Single precision complex complementary error function.                    */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Erfc(rssringoccs_ComplexFloat z)
{
    /*  Declare all necessary variales. C89 requires this at the top.         */
    float x, y;
    rssringoccs_ComplexDouble z_double, erfc_z_double;
    rssringoccs_ComplexFloat erfc_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_CFloat_Real_Part(z);
    y = rssringoccs_CFloat_Imag_Part(z);

    /*  Convert this to a double version.                                     */
    z_double = rssringoccs_CDouble_Rect((double)x, (double)y);

    /*  Compute Erfc(z) on the double version.                                */
    erfc_z_double = rssringoccs_CDouble_Erfc(z_double);

    /*  Extract the real and imaginary parts from erf_z_double, converting    *
     *  them to single precision floats.                                      */
    x = (float)rssringoccs_CDouble_Real_Part(erfc_z_double);
    y = (float)rssringoccs_CDouble_Imag_Part(erfc_z_double);

    /*  Create the float version of erf_z_double and return.                  */
    erfc_z = rssringoccs_CFloat_Rect(x, y);
    return erfc_z;
}
/*  End of rssringoccs_CFloat_Erfc.                                           */

/*  Long double precision complex complementary error function.               */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Erfc(rssringoccs_ComplexLongDouble z)
{
    /*  Declare all necessary variales. C89 requires this at the top.         */
    long double x, y;
    rssringoccs_ComplexDouble z_double, erfc_z_double;
    rssringoccs_ComplexLongDouble erfc_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_CLDouble_Real_Part(z);
    y = rssringoccs_CLDouble_Imag_Part(z);

    /*  Convert this to a double version.                                     */
    z_double = rssringoccs_CDouble_Rect((double)x, (double)y);

    /*  Compute Erfc(z) on the double version.                                */
    erfc_z_double = rssringoccs_CDouble_Erfc(z_double);

    /*  Extract the real and imaginary parts from erf_z_double, converting    *
     *  them to long double precision.                                        */
    x = (long double)rssringoccs_CDouble_Real_Part(erfc_z_double);
    y = (long double)rssringoccs_CDouble_Imag_Part(erfc_z_double);

    /*  Create the long double version of erf_z_double and return.            */
    erfc_z = rssringoccs_CLDouble_Rect(x, y);
    return erfc_z;
}
/*  End of rssringoccs_CFloat_Erfc.                                           */
