/******************************************************************************
 *                                Bessel J0                                   *
 ******************************************************************************
 *  This file contains functions for computing the Bessel J0 function.        * 
 ******************************************************************************
 * We define J_0(x) as the power series solution to the ODE:                  *
 *                                                                            *
 *      x^2y''(x) + xy'y(x) + x^2y(x) = 0                                     *
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
 *                      ___                                                 *
 *                      \      (-1)^n x^2n /                                  *
 *         J_0(x)  =    /__             / (n)!^2 * 4^n                      *
 *                      n = 0                                                 *
 *                                                                            *
 *      For large arguments the asymptotic expansion is used. This is defined *
 *      by the following series:                                              *
 *                                                                            *
 *                      ___                                                 *
 *                      \      cos(z) a_{2n} /    + sin(z) a_{2n+1} /         *
 *          J_0(x)  ~   /__               / x^2n                 / x^{2n+1} *
 *                      n = 0                                                 *
 *                                                                            *
 *      Where:                                                                *
 *          a_n = (2n)!^2 /                                                   *
 *                       / 32^n (k!)^3                                        *
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

/*  Coefficients for the Taylor series and asymptotic expansions found here.  */
#include "_math_constants.h"

/*  Main header for the math functions. Contains <math.h> as well.            */
#include "special_functions.h"

/*  Compute the Bessel J_0 function for a floating point number x.            */
float BesselJ0_Float(float x)
{
    /*  Declare necessary variables.                                          */
    float arg = x*x;

    /*  For small arguments, use the Taylor series of J_0.                    */
    if (arg < 50.0){
        float bessel_J0;
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
        return bessel_J0;
    }
    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32) {

        /*  Declare variables used in the asymptotic expansion.               */
        float sinarg, cosarg;

        /*  J_0 is an even function so use the absolute value of x.           */
        x = fabsf(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= sinf(x - PI_BY_FOUR)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_06 + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= cosf(x - PI_BY_FOUR);

        /*  Multiply the result by the coefficient and return.                */
        return (cosarg + sinarg)*SQRT_2_BY_PI/sqrtf(x);

    }
    /*  For very large arguments, use the limit (which is zero).              */
    else return 0.0;
}

/*  Compute the Bessel J_0 function for a double precision number x.          */
double BesselJ0_Double(double x)
{
    /*  Declare necessary variables.                                          */
    double arg = x*x;

    /*  For small arguments, use the Taylor series of J_0.                    */
    if (arg < 150.0){
        double bessel_J0;
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
        return bessel_J0;
    }
    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32) {

        /*  Declare variables used in the asymptotic expansion.               */
        double sinarg, cosarg;

        /*  J_0 is an even function so use the absolute value of x.           */
        x = fabs(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= sin(x - PI_BY_FOUR)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_06 + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= cos(x - PI_BY_FOUR);

        /*  Multiply the result by the coefficient and return.                */
        return (cosarg + sinarg)*SQRT_2_BY_PI/sqrt(x);

    }
    /*  For very large arguments, use the limit (which is zero).              */
    else return 0.0;
}

/*  Compute the Bessel I_0 function for a long double precision number x.     */
long double BesselJ0_Long_Double(long double x)
{
    /*  Declare necessary variables.                                          */
    long double arg = x*x;

    /*  For small arguments, use the Taylor series of J_0.                    */
    if (arg < 150.0){
        long double bessel_J0;
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
        return bessel_J0;
    }
    /*  For large arguments use the asymptotic expansion.                     */
    else if (arg < 1.0e32) {

        /*  Declare variables used in the asymptotic expansion.               */
        double sinarg, cosarg;

        /*  J_0 is an even function so use the absolute value of x.           */
        x = fabsl(x);

        /*  The argument for the asymptotic expansion is 1/x^2.               */
        arg = 1.0/arg;

        /*  Use Horner's method to compute the polynomial part.               */
        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;

        /*  Multiply the output by the coefficient factor.                    */
        sinarg *= sinl(x - PI_BY_FOUR)/x;

        /*  Do the same as above for the Cosine portion.                      */
        cosarg  = arg * BESSEL_J0_ASYM_08 + BESSEL_J0_ASYM_06;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= cosl(x - PI_BY_FOUR);

        /*  For very large arguments, use the limit (which is zero).          */
        return (cosarg + sinarg)*SQRT_2_BY_PI/sqrtl(x);

    }
    /*  For very large arguments, use the limit (which is zero).              */
    else return 0.0;
}

/*  For all integer types, convert to double and compute.                     */
RSSRINGOCCSNonFloatInputForFloatOutput(BesselJ0);
