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
 *                                Bessel J0                                   *
 ******************************************************************************
 *  This file contains functions for computing the Bessel J0 function.        *
 ******************************************************************************
 * We define J_0(x) as the power series solution to the ODE:                  *
 *                                                                            *
 *           d^2 y(x)      dy(x)                                              *
 *      x^2 ---------  + x ----- + x^2 y(x) = 0                               *
 *             dx^2         dx                                                *
 *                                                                            *
 *  Where J_0 is the Bessel function of the First kind with alpha = 0.        *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  BesselJ0                                                                  *
 *  Purpose:                                                                  *
 *      Compute the J_0 bessel function for a real argument.                  *
 *  Arguments:                                                                *
 *      x (float, double, or long double):                                    *
 *          A real number, the argument for J_0(x).                           *
 *  Output:                                                                   *
 *      bessel_J0:                                                            *
 *          The Bessel function J_0(x).                                       *
 *  Method:                                                                   *
 *      For small values, the Taylor expansion is used. The J_0(x) function   *
 *      can be defined by the following series:                               *
 *                                                                            *
 *                      ___                                                   *
 *                      \       (-1)^n x^2n                                   *
 *         J_0(x)  =    /__     ------------                                  *
 *                     n = 0    (n)!^2 * 4^n                                  *
 *                                                                            *
 *      For large arguments the asymptotic expansion is used. This is defined *
 *      by the following series:                                              *
 *                                                                            *
 *                      ___          _                                  _     *
 *                      \           |  cos(z) a_{2n}    sin(z) a_{2n+1}  |    *
 *          J_0(x)  ~   /__  (-1)^n |  ------------- -  ---------------  |    *
 *                     n = 1         -      x^2n            x^{2n+1}    -     *
 *                                                                            *
 *      Where:                                                                *
 *                         (2n)!^2                                            *
 *          a_n = (-1)^n ------------                                         *
 *                        32^n (n!)^3                                         *
 *                                                                            *
 *          z   = x - pi/4                                                    *
 *                                                                            *
 *      Note that this expansion diverges for all real numbers. To make use   *
 *      of this series we must stop the sum at a particular N. We compute     *
 *      between 8 and 10 terms of this expansion, depending on the precision  *
 *      desired (float, double, long double).                                 *
 *  Error:                                                                    *
 *      In the region in which the Taylor series is used, relative error is   *
 *      10^-16. In the hand-off region with the asymptotic expansion, the     *
 *      error is 10^-9 but quickly drops back to 10^-16.                      *
 *      The regions where the Taylor series is used are listed below:         *
 *          float:           (-7.07,  7.07)                                   *
 *          double:         (-12.24, 12.24)                                   *
 *          long double:    (-12.24, 12.24)                                   *
 *      The alternating series test gives the error of the partial sums of    *
 *      the Taylor expansion. This, combined with trial and error, produced   *
 *      these selected ranges.                                                *
 ******************************************************************************/

/*  The C standard library header math.h is included here, as are aliases for *
 *  various functions, the macros INFINITY and NAN, as well as the max legal  *
 *  values for the exponential function which don't return INFINITY.          */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Complex variables and functions defined here.                             */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Taylor Expansion of Bessel J_0(x).                                        */
static float taylorf[17] = {
     1.0F,
    -0.25F,
     1.56250e-2F,
    -4.340277777F,
     6.781684027F,
    -6.781684027F,
     4.709502797F,
    -2.402807549F,
     9.385966990F,
    -2.896903392F,
     7.242258480F,
    -1.496334396F,
     2.597802772F,
    -3.842903509F,
     4.901662639F,
    -5.446291821F,
     5.318644356F
};

static double taylor[31] = {
     1.0,
    -0.25,
     1.56250e-2,
    -4.34027777777777777e-4,
     6.78168402777777777e-6,
    -6.78168402777777777e-8,
     4.70950279706790123e-10,
    -2.40280754952443940e-12,
     9.38596699032984142e-15,
    -2.89690339207711155e-17,
     7.24225848019277887e-20,
    -1.49633439673404522e-22,
     2.59780277210771740e-25,
    -3.84290350903508491e-28,
     4.90166263907536340e-31,
    -5.44629182119484823e-34,
     5.31864435663559397e-37,
    -4.60090342269515049e-40,
     3.55007980146230748e-43,
    -2.45850401763317692e-46,
     1.53656501102073557e-49,
    -8.71068600351890918e-53,
     4.49932128280935391e-56,
    -2.12633330945621640e-59,
     9.22887721118149482e-63,
    -3.69155088447259792e-66,
     1.36521852236412645e-69,
    -4.68181934967121554e-73,
     1.49292708854311720e-76,
    -4.43795210625183472e-80,
     1.23276447395884297e-83
};


static long double taylorl[31] = {
     1.0L,
    -0.25L,
     1.56250e-2L,
    -4.34027777777777777777777777778e-4L,
     6.78168402777777777777777777778e-6L,
    -6.78168402777777777777777777778e-8L,
     4.70950279706790123456790123457e-10L,
    -2.40280754952443940539178634417e-12L,
     9.38596699032984142731166540690e-15L,
    -2.89690339207711155163940290337e-17L,
     7.24225848019277887909850725841e-20L,
    -1.49633439673404522295423703686e-22L,
     2.59780277210771740096221707789e-25L,
    -3.84290350903508491266600159451e-28L,
     4.90166263907536340901275713585e-31L,
    -5.44629182119484823223639681761e-34L,
     5.31864435663559397679335626720e-37L,
    -4.60090342269515049895619054256e-40L,
     3.55007980146230748376249270259e-43L,
    -2.45850401763317692781336059736e-46L,
     1.53656501102073557988335037335e-49L,
    -8.71068600351890918301219032512e-53L,
     4.49932128280935391684513963074e-56L,
    -2.12633330945621640682662553438e-59L,
     9.22887721118149482129611777074e-63L,
    -3.69155088447259792851844710830e-66L,
     1.36521852236412645285445529153e-69L,
    -4.68181934967121554476836519729e-73L,
     1.49292708854311720177562665730e-76L,
    -4.43795210625183472584906854132e-80L,
     1.23276447395884297940251903925e-83L
};

/*  Asympotic Expansion of Bessel J_0(x).                                     */
#define BESSEL_J0_ASYM_00_F  1.0F
#define BESSEL_J0_ASYM_01_F  0.1250F
#define BESSEL_J0_ASYM_02_F -0.07031250F
#define BESSEL_J0_ASYM_03_F -0.07324218750F
#define BESSEL_J0_ASYM_04_F  0.1121520996093750F
#define BESSEL_J0_ASYM_05_F  0.2271080017089843750F
#define BESSEL_J0_ASYM_06_F -0.57250142097473144531250F
#define BESSEL_J0_ASYM_07_F -1.72772750258445739746093750F
#define BESSEL_J0_ASYM_08_F  6.07404200127348303794860839844F
#define BESSEL_J0_ASYM_09_F  24.3805296995560638606548309326F
#define BESSEL_J0_ASYM_10_F -110.017140269246738171204924583F
#define BESSEL_J0_ASYM_11_F -551.335896122020585607970133424F

#define BESSEL_J0_ASYM_00  1.0
#define BESSEL_J0_ASYM_01  0.1250
#define BESSEL_J0_ASYM_02 -0.07031250
#define BESSEL_J0_ASYM_03 -0.07324218750
#define BESSEL_J0_ASYM_04  0.1121520996093750
#define BESSEL_J0_ASYM_05  0.2271080017089843750
#define BESSEL_J0_ASYM_06 -0.57250142097473144531250
#define BESSEL_J0_ASYM_07 -1.72772750258445739746093750
#define BESSEL_J0_ASYM_08  6.07404200127348303794860839844
#define BESSEL_J0_ASYM_09  24.3805296995560638606548309326
#define BESSEL_J0_ASYM_10 -110.017140269246738171204924583
#define BESSEL_J0_ASYM_11 -551.335896122020585607970133424

#define BESSEL_J0_ASYM_00_L  1.0L
#define BESSEL_J0_ASYM_01_L  0.1250L
#define BESSEL_J0_ASYM_02_L -0.07031250L
#define BESSEL_J0_ASYM_03_L -0.07324218750L
#define BESSEL_J0_ASYM_04_L  0.1121520996093750L
#define BESSEL_J0_ASYM_05_L  0.2271080017089843750L
#define BESSEL_J0_ASYM_06_L -0.57250142097473144531250L
#define BESSEL_J0_ASYM_07_L -1.72772750258445739746093750L
#define BESSEL_J0_ASYM_08_L  6.07404200127348303794860839844L
#define BESSEL_J0_ASYM_09_L  24.3805296995560638606548309326L
#define BESSEL_J0_ASYM_10_L -110.017140269246738171204924583L
#define BESSEL_J0_ASYM_11_L -551.335896122020585607970133424L


/*  Compute the Bessel J_0 function for a floating point number x.            */
float rssringoccs_Float_Bessel_J0(float x)
{
    /*  Declare necessary variables. C89 requires declaring these at the top. */
    float bessel_J0, arg;
    float sinarg, cosarg;

    /*  Bessel J0 is even and in terms of the square of x, so compute this.   */
    arg = x*x;

    /*  For small arguments, use the Taylor series of J_0.                    */
    if (arg < 4.0F)
        bessel_J0 = rssringoccs_Real_Poly_Float_Coeffs(taylorf, 7U, arg);
    else if (arg < 16.0F)
        bessel_J0 = rssringoccs_Real_Poly_Float_Coeffs(taylorf, 10U, arg);
    else if (arg < 25.0F)
        bessel_J0 = rssringoccs_Real_Poly_Float_Coeffs(taylorf, 12U, arg);
    else if (arg < 36.0F)
        bessel_J0 = rssringoccs_Real_Poly_Float_Coeffs(taylorf, 13U, arg);
    else if (arg < 49.0F)
        bessel_J0 = rssringoccs_Real_Poly_Float_Coeffs(taylorf, 15U, arg);
    else if (arg < 64.0F)
        bessel_J0 = rssringoccs_Real_Poly_Float_Coeffs(taylorf, 16U, arg);

    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32F)
    {
        /*  J_0 is an even function so use the absolute value of x.           */
        x = rssringoccs_Float_Abs(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0F/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07_F + BESSEL_J0_ASYM_05_F;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03_F;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01_F;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= rssringoccs_Float_Sin(x - rssringoccs_Pi_By_Four_F)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_06_F + BESSEL_J0_ASYM_04_F;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02_F;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00_F;
        cosarg *= rssringoccs_Float_Cos(x - rssringoccs_Pi_By_Four_F);

        /*  Multiply the result by the coefficient and return.                */
        bessel_J0 = (cosarg + sinarg)*rssringoccs_Sqrt_Two_By_Pi_F;
        bessel_J0 = bessel_J0 / rssringoccs_Float_Sqrt(x);
    }

    /*  For very large arguments, use the limit (which is zero).              */
    else
        bessel_J0 = 0.0F;

    return bessel_J0;
}

/*  Compute the Bessel J_0 function for a double precision number x.          */
double rssringoccs_Double_Bessel_J0(double x)
{
    /*  Declare necessary variables. C89 requires declaring these at the top. */
    double bessel_J0, arg;
    double sinarg, cosarg;

    /*  Bessel J0 is even and in terms of the square of x, so compute this.   */
    arg = x*x;

    /*  For small arguments, use the Taylor series of J_0.                    */
    if (arg < 4.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 12U, arg);
    else if (arg < 16.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 16U, arg);
    else if (arg < 25.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 18U, arg);
    else if (arg < 36.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 19U, arg);
    else if (arg < 49.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 21U, arg);
    else if (arg < 64.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 23U, arg);
    else if (arg < 81.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 24U, arg);
    else if (arg < 100.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 26U, arg);
    else if (arg < 121.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 27U, arg);
    else if (arg < 144.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 29U, arg);
    else if (arg < 196.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 30U, arg);
    else if (arg < 225.0)
        bessel_J0 = rssringoccs_Real_Poly_Double_Coeffs(taylor, 31U, arg);

    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32)
    {
        /*  J_0 is an even function so use the absolute value of x.           */
        x = rssringoccs_Double_Abs(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= rssringoccs_Double_Sin(x - rssringoccs_Pi_By_Four)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_06 + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= rssringoccs_Double_Cos(x - rssringoccs_Pi_By_Four);

        /*  Multiply the result by the coefficient and return.                */
        bessel_J0 = (cosarg + sinarg)*rssringoccs_Sqrt_Two_By_Pi;
        bessel_J0 = bessel_J0 / rssringoccs_Double_Sqrt(x);
    }

    /*  For very large arguments, use the limit (which is zero).              */
    else
        bessel_J0 = 0.0;

    return bessel_J0;
}

/*  Compute the Bessel I_0 function for a long double precision number x.     */
long double rssringoccs_LDouble_Bessel_J0(long double x)
{
    /*  Declare necessary variables. C89 requires declaring these at the top. */
    long double bessel_J0, arg;
    long double sinarg, cosarg;

    /*  Bessel J0 is even and in terms of the square of x, so compute this.   */
    arg = x*x;

    /*  For small arguments, use the Taylor series of J_0.                    */
    if (arg < 4.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 12U, arg);
    else if (arg < 16.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 16U, arg);
    else if (arg < 25.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 18U, arg);
    else if (arg < 36.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 19U, arg);
    else if (arg < 49.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 21U, arg);
    else if (arg < 64.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 23U, arg);
    else if (arg < 81.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 24U, arg);
    else if (arg < 100.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 26U, arg);
    else if (arg < 121.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 27U, arg);
    else if (arg < 144.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 29U, arg);
    else if (arg < 169.0L)
        bessel_J0 = rssringoccs_Real_Poly_LDouble_Coeffs(taylorl, 30U, arg);

    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32L)
    {
        /*  J_0 is an even function so use the absolute value of x.           */
        x = rssringoccs_LDouble_Abs(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0L/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07_L + BESSEL_J0_ASYM_05_L;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03_L;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01_L;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= rssringoccs_LDouble_Sin(x - rssringoccs_Pi_By_Four_L)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_08_L + BESSEL_J0_ASYM_06_L;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_04_L;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02_L;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00_L;
        cosarg *= rssringoccs_LDouble_Cos(x - rssringoccs_Pi_By_Four_L);

        /*  For very large arguments, use the limit (which is zero).          */
        bessel_J0 = (cosarg + sinarg)*rssringoccs_Sqrt_Two_By_Pi_L;
        bessel_J0 = bessel_J0 / rssringoccs_LDouble_Sqrt(x);
    }

    /*  For very large arguments, use the limit (which is zero).              */
    else
        bessel_J0 = 0.0L;

    return bessel_J0;
}
