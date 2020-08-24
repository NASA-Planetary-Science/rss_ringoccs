/******************************************************************************
 *                                Bessel I0                                   *
 ******************************************************************************
 *  This file contains functions for computing the Bessel I0 function.        * 
 ******************************************************************************
 * We define I_0(x) as follows:                                               *
 *                                                                            *
 *      I_0(x)  =  J_0(ix)                                                    *
 *                                                                            *
 *  Where J_0 is the Bessel function of the First kind with alpha = 0.        *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  BesselI0:                                                                 *
 *  Purpose:                                                                  *
 *      Compute the I_0 bessel function for a real argument.                  *
 *  Arguments:                                                                *
 *      x (real number, float, int, char, double, etc.):                      *
 *          A real number, the argument for J_0(x).                           *
 *  Output:                                                                   *
 *      bessel_I0:                                                            *
 *          The Bessel function I_0(x).                                       *
 *  Method:                                                                   *
 *      This uses the standard Taylor series for small inputs, and            *
 *      asymptotic expansions for large inputs. The Taylor series is:         *
 *                                                                            *
 *                   ___                                                    *
 *                   \      x^(2n)/                                           *
 *        I0(x)=     /__       / ((2n)! * 2^n)^2                            *
 *                   n = 0                                                    *
 *                                                                            *
 *      This can be obtained by simply evaluating ix into the J_0 function.   *
 *                                                                            *
 *      The asymptotic expansions can be obtained by iteratively using        *
 *      integration by parts. This series diverge for all x, but by cutting   *
 *      off the expansion at a particular N, one finds a very good            *
 *      approximation. For large values, we have:                             *
 *                                                                            *
 *          I_0(x)  ~   exp(x) / sqrt(2 * pi * x)                             *
 *                                                                            *
 *      For very large values, we simply return infinity.                     *
 *  Error:                                                                    *
 *      In the region in which the Taylor series is used, relative error is   *
 *      10^-16. In the hand-off region with the asymptotic expansion, the     *
 *      error is 10^-9 but quickly drops back to 10^-16.                      *
 *      The regions where the Taylor series is used are listed below:         *
 *          float:           (-3.46, 3.46)                                    *
 *          double:          (-4.00, 4.00)                                    *
 *          long double:     (-4.35, 4.35)                                    *
 *  Notes:                                                                    *
 *      This code was written with an emphasis on clarity and accessibility   *
 *      to a wide audience without the need for anything too sophisticated.   *
 *      More accurate methods involving rational functions and Chebyshev      *
 *      polynomials exist. See the GNU Scientific Library for source code.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 15, 2019                                            *
 ******************************************************************************/

/*  Coefficients for the Taylor series and asymptotic expansions found here.  */
#include "_math_constants.h"

/*  Main header for the math functions. Contains <math.h> as well.            */
#include "special_functions.h"

/*  Compute the Bessel I_0 function for a floating point number x.            */
float BesselI0_Float(float x)
{
    /*  Declare necessary variables.                                          */
    float bessel_I0, arg;

    /*  I_0 is an even function, so compute the absolute value of x.          */
    x = fabsf(x);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (x < 12.0){

        /*  The series is in powers of x^2, so use Horner's Method with that. */
        arg = x*x;
        bessel_I0 = arg * BESSEL_I0_TAYLOR_16 + BESSEL_I0_TAYLOR_15;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_14;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_13;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_12;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_11;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_10;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_09;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_08;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_07;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_06;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_05;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_04;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_00;
        return bessel_I0;
    }
    /*  For larger values, use the asymptotic expansion.                      */
    else if (x < 87.49) {

        /*  The asymptotic expansion is in terms of 1/x.                      */
        arg = 1.0/x;

        /*  Compute the polynomial term using Horner's Method.                */
        bessel_I0 = arg * BESSEL_I0_ASYM_04 + BESSEL_I0_ASYM_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_00;

        /*  Multiply by the coefficient factor and return.                    */
        return bessel_I0 * expf(x)/sqrtf(TWO_PI*x);
    }
    /*  For very large inputs, return INFINITY (standard in C99).             */
    else return INFINITY;
}

/*  Compute the Bessel I_0 function for a double precision number x.          */
double BesselI0_Double(double x)
{
    /*  Declare necessary variables.                                          */
    double bessel_I0, arg;

    /*  I_0 is symmetric so compute the absolute value of x and use that.     */
    x = fabs(x);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (x < 16.0){

        /*  The series is in powers of x^2, so use Horner's Method with that. */
        arg = x*x;
        bessel_I0 = arg * BESSEL_I0_TAYLOR_22 + BESSEL_I0_TAYLOR_21;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_20;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_19;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_18;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_17;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_16;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_15;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_14;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_13;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_12;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_11;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_10;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_09;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_08;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_07;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_06;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_05;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_04;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_00;
        return bessel_I0;
    }
    /*  For larger values, use the asymptotic expansion.                      */
    else if (x < 709.0) {

        /*  The asymptotic expansion is in terms of 1/x.                      */
        arg = 1.0/x;

        /*  Compute the polynomial term using Horner's Method.                */
        bessel_I0 = arg * BESSEL_I0_ASYM_06 + BESSEL_I0_ASYM_05;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_04;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_00;

        /*  Multiply by the coefficient factor and return.                    */
        return bessel_I0 * exp(x)/sqrt(TWO_PI*x);

    }
    /*  For very large inputs, return INFINITY (standard in C99).             */
    else return INFINITY;
}

/*  Compute the Bessel I_0 function for a long double precision number x.     */
long double BesselI0_Long_Double(long double x)
{
    /*  Declare necessary variables.                                          */
    long double bessel_I0, arg;

    /*  I_0 is symmetric so compute the absolute value of x and use that.     */
    x = fabsl(x);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (x < 19.0){

        /*  The series is in powers of x^2, so use Horner's Method with that. */
        arg = x*x;
        bessel_I0 = arg * BESSEL_I0_TAYLOR_24 + BESSEL_I0_TAYLOR_23;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_22;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_21;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_20;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_19;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_18;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_17;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_16;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_15;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_14;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_13;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_12;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_11;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_10;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_09;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_08;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_07;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_06;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_05;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_04;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_TAYLOR_00;
        return bessel_I0;
    }
    /*  For larger values, use the asymptotic expansion.                      */
    else if (x < 11356.34) {

        /*  The asymptotic expansion is in terms of 1/x.                      */
        arg = 1.0/x;

        /*  Compute the polynomial term using Horner's Method.                */
        bessel_I0 = arg * BESSEL_I0_ASYM_06 + BESSEL_I0_ASYM_05;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_04;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_00;

        /*  Multiply by the coefficient factor and return.                    */
        return bessel_I0 * expl(x)/sqrtl(TWO_PI*x);

    }
    /*  For very large inputs, return INFINITY (standard in C99).             */
    else return INFINITY;
}

/*  For all integer types, convert to double and compute.                     */
RSSRINGOCCSNonFloatInputForFloatOutput(BesselI0);
