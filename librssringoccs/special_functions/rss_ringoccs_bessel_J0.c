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
#define BESSEL_J0_TAYLOR_00  1.0
#define BESSEL_J0_TAYLOR_01 -0.25
#define BESSEL_J0_TAYLOR_02  1.56250e-2
#define BESSEL_J0_TAYLOR_03 -4.34027777777777777777777777778e-4
#define BESSEL_J0_TAYLOR_04  6.78168402777777777777777777778e-6
#define BESSEL_J0_TAYLOR_05 -6.78168402777777777777777777778e-8
#define BESSEL_J0_TAYLOR_06  4.70950279706790123456790123457e-10
#define BESSEL_J0_TAYLOR_07 -2.40280754952443940539178634417e-12
#define BESSEL_J0_TAYLOR_08  9.38596699032984142731166540690e-15
#define BESSEL_J0_TAYLOR_09 -2.89690339207711155163940290337e-17
#define BESSEL_J0_TAYLOR_10  7.24225848019277887909850725841e-20
#define BESSEL_J0_TAYLOR_11 -1.49633439673404522295423703686e-22
#define BESSEL_J0_TAYLOR_12  2.59780277210771740096221707789e-25
#define BESSEL_J0_TAYLOR_13 -3.84290350903508491266600159451e-28
#define BESSEL_J0_TAYLOR_14  4.90166263907536340901275713585e-31
#define BESSEL_J0_TAYLOR_15 -5.44629182119484823223639681761e-34
#define BESSEL_J0_TAYLOR_16  5.31864435663559397679335626720e-37
#define BESSEL_J0_TAYLOR_17 -4.60090342269515049895619054256e-40
#define BESSEL_J0_TAYLOR_18  3.55007980146230748376249270259e-43
#define BESSEL_J0_TAYLOR_19 -2.45850401763317692781336059736e-46
#define BESSEL_J0_TAYLOR_20  1.53656501102073557988335037335e-49
#define BESSEL_J0_TAYLOR_21 -8.71068600351890918301219032512e-53
#define BESSEL_J0_TAYLOR_22  4.49932128280935391684513963074e-56
#define BESSEL_J0_TAYLOR_23 -2.12633330945621640682662553438e-59
#define BESSEL_J0_TAYLOR_24  9.22887721118149482129611777074e-63
#define BESSEL_J0_TAYLOR_25 -3.69155088447259792851844710830e-66
#define BESSEL_J0_TAYLOR_26  1.36521852236412645285445529153e-69
#define BESSEL_J0_TAYLOR_27 -4.68181934967121554476836519729e-73
#define BESSEL_J0_TAYLOR_28  1.49292708854311720177562665730e-76
#define BESSEL_J0_TAYLOR_29 -4.43795210625183472584906854132e-80
#define BESSEL_J0_TAYLOR_30  1.23276447395884297940251903925e-83

/*  Asympotic Expansion of Bessel J_0(x).                                     */
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


/*  Compute the Bessel J_0 function for a floating point number x.            */
float rssringoccs_Float_Bessel_J0(float x)
{
    /*  Declare necessary variables. C89 requires declaring these at the top. */
    float bessel_J0, arg;
    float sinarg, cosarg;

    /*  Bessel J0 is even and in terms of the square of x, so compute this.   */
    arg = x*x;

    /*  For small arguments, use the Taylor series of J_0.                    */
    if (arg < 50.0)
    {
        bessel_J0 = arg * BESSEL_J0_TAYLOR_15 + BESSEL_J0_TAYLOR_14;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_13;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_12;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_11;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_10;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_09;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_08;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_07;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_06;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_05;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_04;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_03;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_02;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_01;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_00;
    }

    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32)
    {
        /*  J_0 is an even function so use the absolute value of x.           */
        x = rssringoccs_Float_Abs(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= rssringoccs_Float_Sin(x - rssringoccs_Pi_By_Four)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_06 + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= rssringoccs_Float_Cos(x - rssringoccs_Pi_By_Four);

        /*  Multiply the result by the coefficient and return.                */
        bessel_J0 = (cosarg + sinarg)*rssringoccs_Sqrt_Two_By_Pi;
        bessel_J0 = bessel_J0 / rssringoccs_Float_Sqrt(x);
    }

    /*  For very large arguments, use the limit (which is zero).              */
    else
        bessel_J0 = 0.0;

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
    if (arg < 150.0)
    {
        bessel_J0 = arg * BESSEL_J0_TAYLOR_22 + BESSEL_J0_TAYLOR_21;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_20;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_19;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_18;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_17;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_16;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_15;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_14;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_13;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_12;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_11;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_10;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_09;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_08;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_07;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_06;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_05;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_04;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_03;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_02;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_01;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_00;
    }

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
    if (arg < 150.0)
    {
        bessel_J0 = arg * BESSEL_J0_TAYLOR_24 + BESSEL_J0_TAYLOR_23;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_22;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_21;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_20;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_19;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_18;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_17;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_16;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_15;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_14;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_13;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_12;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_11;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_10;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_09;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_08;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_07;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_06;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_05;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_04;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_03;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_02;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_01;
        bessel_J0 = arg * bessel_J0 + BESSEL_J0_TAYLOR_00;
    }

    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32)
    {
        /*  J_0 is an even function so use the absolute value of x.           */
        x = rssringoccs_LDouble_Abs(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= rssringoccs_LDouble_Sin(x - rssringoccs_Pi_By_Four)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_08 + BESSEL_J0_ASYM_06;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= rssringoccs_LDouble_Cos(x - rssringoccs_Pi_By_Four);

        /*  For very large arguments, use the limit (which is zero).          */
        bessel_J0 = (cosarg + sinarg)*rssringoccs_Sqrt_Two_By_Pi;
        bessel_J0 = bessel_J0 / rssringoccs_LDouble_Sqrt(x);
    }

    /*  For very large arguments, use the limit (which is zero).              */
    else
        bessel_J0 = 0.0;

    return bessel_J0;
}
