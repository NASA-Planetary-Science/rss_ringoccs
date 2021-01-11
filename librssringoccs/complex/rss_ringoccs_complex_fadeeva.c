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
 *                     rss_ringoccs_complex_faddeeva                          *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the Faddeeva function.                   *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CFloat_Faddeeva:                                          *
 *      rssringoccs_CDouble_Faddeeva:                                         *
 *      rssringoccs_CLDouble_Faddeeva:                                        *
 *  Purpose:                                                                  *
 *      Computes the Faddeeva function w(z):                                  *
 *                                                                            *
 *          w(z) = exp(-z^2) Erfc(-iz)                                        *
 *                                                                            *
 *      Where Erfc is the complementary error function.                       *
 *  Arguments:                                                                *
 *      z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      w (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          The Faddeeva function w(z).                                       *
 *  Method:                                                                   *
 *      Analyze the input to see which region of the complex plane it lies in *
 *      and use the fadeeva and Erfcx functions to compute.                   *
 *  NOTES:                                                                    *
 *      No actual float or long double algorithms have been implemented by    *
 *      rss_ringoccs. We simply cast to double and compute using that version.*
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
 *  Date:       December 28, 2020                                             *
 ******************************************************************************/

/******************************************************************************
 *                             LIBCERF LICENSE                                *
 ******************************************************************************/

/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File w_of_z.c:
 *   Computation of Faddeeva's complex scaled error function,
 *      w(z) = exp(-z^2) * erfc(-i*z),
 *   nameless function (7.1.3) of Abramowitz&Stegun (1964),
 *   also known as the plasma dispersion function.
 *
 *   This implementation uses a combination of different algorithms.
 *   See man 3 w_of_z for references.
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
 * Man page:
 *   w_of_z(3)
 */


#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/******************************************************************************
 *  Precomputed table of expa2n2[n-1] = exp(-a2*n*n)                          *
 *  for double-precision a2 = 0.26865... in w_of_z, below.                    *
 ******************************************************************************/

static const double expa2n2[] = {
    7.64405281671221563e-01,
    3.41424527166548425e-01,
    8.91072646929412548e-02,
    1.35887299055460086e-02,
    1.21085455253437481e-03,
    6.30452613933449404e-05,
    1.91805156577114683e-06,
    3.40969447714832381e-08,
    3.54175089099469393e-10,
    2.14965079583260682e-12,
    7.62368911833724354e-15,
    1.57982797110681093e-17,
    1.91294189103582677e-20,
    1.35344656764205340e-23,
    5.59535712428588720e-27,
    1.35164257972401769e-30,
    1.90784582843501167e-34,
    1.57351920291442930e-38,
    7.58312432328032845e-43,
    2.13536275438697082e-47,
    3.51352063787195769e-52,
    3.37800830266396920e-57,
    1.89769439468301000e-62,
    6.22929926072668851e-68,
    1.19481172006938722e-73,
    1.33908181133005953e-79,
    8.76924303483223939e-86,
    3.35555576166254986e-92,
    7.50264110688173024e-99,
    9.80192200745410268e-106,
    7.48265412822268959e-113,
    3.33770122566809425e-120,
    8.69934598159861140e-128,
    1.32486951484088852e-135,
    1.17898144201315253e-143,
    6.13039120236180012e-152,
    1.86258785950822098e-160,
    3.30668408201432783e-169,
    3.43017280887946235e-178,
    2.07915397775808219e-187,
    7.36384545323984966e-197,
    1.52394760394085741e-206,
    1.84281935046532100e-216,
    1.30209553802992923e-226,
    5.37588903521080531e-237,
    1.29689584599763145e-247,
    1.82813078022866562e-258,
    1.50576355348684241e-269,
    7.24692320799294194e-281,
    2.03797051314726829e-292,
    3.34880215927873807e-304,

    /*  Underflow (also prevents reads past array end, below).                */
    0.0
};

/*  Double precsision complex Faddeeva function.                              */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Faddeeva(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. We'll declare all variables used in any  *
     *  branch here for clarity. C89/C90 requires variables to be delcared at *
     *  the top of a code block and does not allow mixed declarations/code.   */
    int n, dn;
    double z_x, z_y, abs_x, abs_y;
    double w_x, w_y;
    double sum1, sum2, sum3, sum4, sum5;
    double xs, abs_y_by_xs, xs_by_abs_y;
    double dr, di;
    double denom, coef;
    double c0, c1, c2, c3, c4;
    double nu, wi, wr;
    double prod2ax, prodm2ax, expx2;
    double sinxy, sin2xy, cos2xy, coef1, coef2, sincxy, sinc2xy;
    double x2, ax2, exp2ax, expm2ax, expx2erfcxy;
    double re_temp, n0, dx;
    double np, nm, tp, tm;
    double exp1, exp1dn;
    rssringoccs_ComplexDouble w, temp;

    /*  Define a to be pi / sqrt(-ln(DBL_EPSILON*0.5)). This value assumes    *
     *  DBL_EPSILON = 2.22045e-16. This is in compliance with IEEE 754 double *
     *  format which uses 52 bits for the mantissa. This yields 2^-52  which  *
     *  is roughly 2.2e-16.                                                   */
    const double a = 0.518321480430085929872;

    /*  Set c to 2*a/pi.                                                      */
    const double c = 0.329973702884629072537;

    /*  And store a^2 as a variable.                                          */
    const double a2 = 0.268657157075235951582;

    /*  Extract the real and imaginary parts from the input z.                */
    z_x = rssringoccs_CDouble_Real_Part(z);
    z_y = rssringoccs_CDouble_Imag_Part(z);

    /*  If we have a purely imaginary input, we get a purely real output. We  *
     *  have:                                                                 *
     *      w(iy) = exp(-(iy)^2) Erfc(-i(iy))                                 *
     *            = exp(y^2) Erfc(y)                                          *
     *  And Erfc of a real number is again real. exp(y^2) Erfc(y) is the      *
     *  definition of the Erfcx function, so use this.                        */
    if (z_x == 0.0)
    {
        w_x = rssringoccs_Double_Erfcx(z_y);
        w_y = 0.0;
        w = rssringoccs_CDouble_Rect(w_x, w_y);
        return w;
    }

    /*  For a purely real input, we have:                                     *
     *      w(x) = exp(-x^2) Erfc(-ix)                                        *
     *  So we need the complementary function for a purely imaginary input.   *
     *  We can use the imaginary part of the Faddeeva function for this.      *
     *  rssringoccs_Double_Faddeeva_Im is defined in rss_ringoccs_math.h.     */
    else if (z_y == 0)
    {
        w_x = rssringoccs_Double_Exp(-z_x*z_x);
        w_y = rssringoccs_Double_Faddeeva_Im(z_x);
        w = rssringoccs_CDouble_Rect(w_x, w_y);
        return w;
    }

    /*  Get the absolute values of both the real and imaginary parts of z.    */
    abs_x = rssringoccs_Double_Abs(z_x);
    abs_y = rssringoccs_Double_Abs(z_y);

    /*  Initialize the output to 0+i0.                                        */
    w = rssringoccs_CDouble_Zero;

    /*  And initialize all of the sum variables to zero.                      */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;

    /*  Original comment from libcerf:                                        *
     *      The continued fraction algorithm is faster. As pointed out by M.  *
     *      Zaghloul the continued fraction seems to give a large relative    *
     *      error in Re w(z) for |x| ~ 6 and small |y|, so use algorithm 816  *
     *      in this region.                                                   */
    if (abs_y > 7.0 || ((abs_x > 6.0) && ((abs_y > 0.1) ||
                                         ((abs_x > 8.0) && (abs_y > 1e-10)) ||
                                         (abs_x > 28.0))))
    {
        /*  Original comment from libcerf:
         *      Poppe & Wijers suggest using a number of terms                *
         *          nu = 3 + 1442 / (26*rho + 77)                             *
         *      where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.   *
         *      (They only use this expansion for rho >= 1, but rho a little  *
         *      less than 1 seems okay too). Instead, I did my own fit to a   *
         *      slightly different function that avoids the hypotenuse        *
         *      calculation, using NLopt to minimize the sum of the squares of*
         *      the errors in nu with the constraint that the estimated nu be *
         *      >= minimum nu to attain machine precision. I also separate the*
         *      regions where nu == 2 and nu == 1.                            */

        /*  Compute for -z if y < 0.                                          */
        if (z_y < 0.0)
            xs = -z_x;
        else
            xs = z_x;

        /*  nu <= 2                                                           */
        if ((abs_x + abs_y) > 4000)
        {
            /*  nu == 1, w(z) = i/sqrt(pi) / z. Scale to avoid overflow.      */
            if ((abs_x + abs_y) > 1e7)
            {
                if (abs_x > abs_y)
                {
                    abs_y_by_xs = abs_y / xs;
                    denom = rssringoccs_Sqrt_One_By_Pi /
                            (xs + abs_y_by_xs*abs_y);
                    w_x = denom*abs_y_by_xs;
                    w_y = denom;
                    w = rssringoccs_CDouble_Rect(w_x, w_y);
                }
                /*  Use the limiting behavior of w(x+iy) if y = inf.          */
                else if (rssringoccs_Is_Inf(abs_y))
                {
                    /*  Check if the x component is a NaN. If it is, or if    *
                     *  if the y component is negative infinity, return       *
                     *  complex NaN.                                          */
                    if ((rssringoccs_Is_NaN(z_x)) || (z_y < 0.0))
                        w = rssringoccs_CDouble_NaN;

                    /*  In the limiting case with y = +infinity we get zero.  */
                    else
                        w = rssringoccs_CDouble_Zero;

                    return w;
                }
                else
                {
                    xs_by_abs_y = xs / abs_y;
                    denom = rssringoccs_Sqrt_One_By_Pi /
                            (xs_by_abs_y*xs + abs_y);
                    w_x = denom;
                    w_y = denom*xs_by_abs_y;
                    w = rssringoccs_CDouble_Rect(w_x, w_y);
                }
            }

            /*  nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5).                 */
            else
            {
                dr = xs*xs - abs_y*abs_y - 0.5;
                di = 2.0*xs*abs_y;
                denom = rssringoccs_Sqrt_One_By_Pi / (dr*dr + di*di);
                w_x = denom * (xs*di - abs_y*dr);
                w_y = denom * (xs*dr + abs_y*di);
                w = rssringoccs_CDouble_Rect(w_x, w_y);
            }
        }

        /*  Compute nu(z) estimate and do general continued fraction.         */
        else
        {
            c0=3.9;
            c1=11.398;
            c2=0.08254;
            c3=0.1421;
            c4=0.2023;
            nu = floor(c0 + c1 / (c2*abs_x + c3*abs_y + c4));
            wr = xs;
            wi = abs_y;
            nu = 0.5 * (nu - 1);
            while (nu > 0.4)
            {
                /*  w <- z - nu/w:                                            */
                denom = nu / (wr*wr + wi*wi);
                wr = xs - wr * denom;
                wi = abs_y + wi * denom;
                nu -= 0.5;
            }
            /*  w(z) = i/sqrt(pi) / w:                                        */
            denom = rssringoccs_Sqrt_One_By_Pi / (wr*wr + wi*wi);
            w_x = denom*wi;
            w_y = denom*wr;
            w = rssringoccs_CDouble_Rect(w_x, w_y);
        }

        if (z_y < 0)
        {
            /*  use w(z) = 2.0*exp(-z*z) - w(-z), but be careful of overflow  *
             *  in exp(-z*z).                                                 */
            w_x = (abs_y - xs) * (xs + abs_y);
            w_y = 2*xs*z_y;
            temp = rssringoccs_CDouble_Rect(w_x, w_y);
            temp = rssringoccs_CDouble_Exp(temp);
            temp = rssringoccs_CDouble_Multiply_Real(2.0, temp);
            w = rssringoccs_CDouble_Subtract(temp, w);

            return w;
        }
        else
            return w;
    }

    /*  Note: The test that seems to be suggested in the paper is             *
     *  x < sqrt(-log(DBL_MIN), about 26.6, since otherwise exp(-x^2)         *
     *  underflows to zero and sum1,sum2,sum4 are zero.  However, long        *
     *  before this occurs, the sum1,sum2,sum4 contributions are              *
     *  negligible in double precision; I find that this happens for x >      *
     *  about 6, for all y.  On the other hand, I find that the case          *
     *  where we compute all of the sums is faster (at least with the         *
     *  precomputed expa2n2 table) until about x=10.  Furthermore, if we      *
     *  try to compute all of the sums for x > 20, I find that we             *
     *  sometimes run into numerical problems because underflow/overflow      *
     *  problems start to appear in the various coefficients of the sums,     *
     *  below.  Therefore, we use x < 10 here.                                */
    else if (abs_x < 10)
    {
        prod2ax = 1.0;
        prodm2ax = 1.0;

        if (rssringoccs_Is_Inf(z_y))
        {
            w = rssringoccs_CDouble_Rect(z_y, z_y);
            return w;
        }

        /*  Compute sum4 and sum5 together as sum5-sum4. This special case is *
         *  needed for accuracy.                                              */
        if (abs_x < 5e-4)
        {
            /*  exp(-x*x) via Taylor compute exp(2*a*x) and exp(-2*a*x) via   *
             *  Taylor, to double precision.                                  */
            x2 = abs_x*abs_x;
            expx2 = 1 - x2 * (1 - 0.5*x2);

            /*  2*a*x.  */
            ax2 = 1.036642960860171859744*abs_x;
            exp2ax = 1 + ax2 * (1 + ax2 * (0.5 + 0.166666666666666666667*ax2));
            expm2ax = 1 - ax2 * (1 - ax2 * (0.5 - 0.166666666666666666667*ax2));

            for (n = 1; ; ++n)
            {
                coef = expa2n2[n-1] * expx2 / (a2*(n*n) + z_y*z_y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum1 += coef;
                sum2 += coef * prodm2ax;
                sum3 += coef * prod2ax;

                /*  really = sum5 - sum4                                      */
                sum5 += coef*(2*a)*n*rssringoccs_Double_Sinh((2*a)*n*abs_x);

                /*  test convergence via sum3.                                */
                if (coef * prod2ax < DBL_EPSILON * sum3)
                    break;
            }
        }

        /*  x > 5e-4, compute sum4 and sum5 separately.                       */
        else
        {
            expx2 = exp(-abs_x*abs_x);
            exp2ax = exp((2*a)*abs_x);
            expm2ax = 1 / exp2ax;
            for (n = 1; ; ++n)
            {
                coef = expa2n2[n-1] * expx2 / (a2*(n*n) + z_y*z_y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum1 += coef;
                sum2 += coef * prodm2ax;
                sum4 += (coef * prodm2ax) * (a*n);
                sum3 += coef * prod2ax;
                sum5 += (coef * prod2ax) * (a*n);
                /*  Test convergence via sum5, since this sum has the slowest *
                 *  decay.                                                    */
                if ((coef * prod2ax) * (a*n) < DBL_EPSILON * sum5)
                    break;
            }
        }

        /*  Avoid spurious overflow for large negative y. For y < -6,         *
         *  erfcx(y) = 2*exp(y*y) to double precision.                        */
        if (z_y > -6.0)
            expx2erfcxy = expx2*rssringoccs_Double_Erfcx(z_y);
        else
            expx2erfcxy = 2*exp(z_y*z_y-abs_x*abs_x);

        /*  Imaginary terms cancel.                                           */
        if (z_y > 5.0)
        {
            sinxy = rssringoccs_Double_Sin(abs_x*z_y);
            sincxy = rssringoccs_Double_Sinc(abs_x*z_y);
            cos2xy = rssringoccs_Double_Cos(2.0*abs_x*z_y);
            w_x = (expx2erfcxy - c*z_y*sum1)*cos2xy +
                  c*abs_x*expx2*sinxy*sincxy;
            w = rssringoccs_CDouble_Rect(w_x, 0.0);
        }
        else
        {
            xs = z_x;
            sinxy = rssringoccs_Double_Sin(xs*z_y);
            sin2xy = rssringoccs_Double_Sin(2.0*xs*z_y);
            cos2xy = rssringoccs_Double_Cos(2.0*xs*z_y);
            coef1 = expx2erfcxy - c*z_y*sum1;
            coef2 = c*xs*expx2;
            sincxy = rssringoccs_Double_Sinc(xs*z_y);
            sinc2xy = rssringoccs_Double_Sinc(2*xs*z_y);
            w_x = coef1 * cos2xy + coef2 * sinxy * sincxy;
            w_y = coef2 * sinc2xy - coef1 * sin2xy;
            w = rssringoccs_CDouble_Rect(w_x, w_y);
        }
    }

    /*  x large: only sum3 & sum5 contribute (see above note).                */
    else
    {
        if (rssringoccs_Is_NaN(abs_x))
        {
            w_x = rssringoccs_NaN;
            w_y = rssringoccs_NaN;
            w = rssringoccs_CDouble_Rect(w_x, w_y);
            return w;
        }
        else if (rssringoccs_Is_NaN(z_y))
        {
            w_x = rssringoccs_NaN;
            w_y = rssringoccs_NaN;
            w = rssringoccs_CDouble_Rect(w_x, w_y);
            return w;
        }

        /*  |y| < 1e-10, so we only need exp(-x*x) term (round instead of     *
         *  ceil as in original paper; note that x/a > 1 here).               */
        re_temp = exp(-z_x*z_x);
        w = rssringoccs_CDouble_Rect(re_temp, 0.0);

        /*  sum in both directions, starting at n0.                           */
        n0 = floor(abs_x/a + 0.5);
        dx = a*n0 - abs_x;
        sum3 = exp(-dx*dx) / (a2*(n0*n0) + z_y*z_y);
        sum5 = a*n0 * sum3;
        exp1 = exp(4*a*dx);
        exp1dn = 1.0;

        /*  Loop over n0-dn and n0+dn terms.                                  */
        for (dn = 1; n0 - dn > 0; ++dn)
        {
            np = n0 + dn;
            nm = n0 - dn;
            tp = exp(-(a*dn+dx)*(a*dn+dx));

            /*  trick to get tm from tp.                                      */
            tm = tp * (exp1dn *= exp1);
            tp /= (a2*(np*np) + z_y*z_y);
            tm /= (a2*(nm*nm) + z_y*z_y);
            sum3 += tp + tm;
            sum5 += a * (np * tp + nm * tm);
            if (a * (np * tp + nm * tm) < DBL_EPSILON * sum5)
                goto finish;
        }
        /*  loop over n0+dn terms only (since n0-dn <= 0).                    */
        while (1)
        {
            np = n0 + dn++;
            tp = exp(-(a*dn+dx)*(a*dn+dx)
            ) / (a2*(np*np) + z_y*z_y);
            sum3 += tp;
            sum5 += a * np * tp;
            if (a * np * tp < DBL_EPSILON * sum5)
                goto finish;
        }
    }

finish:
    {
        w_x = (0.5*c)*z_y*(sum2+sum3);
        w_y = (0.5*c)*rssringoccs_Double_Copysign(sum5-sum4, z_x);
        temp = rssringoccs_CDouble_Rect(w_x, w_y);
        w = rssringoccs_CDouble_Add(w, temp);
        return w;
    }
}
