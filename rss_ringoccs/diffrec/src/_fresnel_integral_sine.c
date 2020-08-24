/*******************************************************************************
 *                               Fresnel Sine                                  *
 *******************************************************************************
 * This C program contains several algorithms for the computation of the       *
 * Fresnel Sine and Cosine integrals. This is meant to test the accuracy and   *
 * efficiency of the various algorithms. This file is not included in the      *
 * setup.py file, and thus will not be compiled. The fresnel sine and cosine   *
 * functions that are used in the _special_functions.c file are the functions  *
 * that had the best combination of speed and accuracy. This file is kept so   *
 * that users may experiment with the various known algorithms.                *
 *******************************************************************************
 * We define the Fresnel Sine and Cosine Integrals as follows:                 *
 *           x                                                                 *
 *           -                                                                 *
 *          | |                                                                *
 * S(x) =   |   sin(t^2) dt                                                    *
 *        | |                                                                  *
 *         -                                                                   *
 *         0                                                                   *
 *******************************************************************************
 * It is very common for a pi/2 to be placed inside the sine and cosine terms, *
 * and thus in translating one would need to scale x by sqrt(2/pi) and scale   *
 * the results by sqrt(pi/2). Several of the algorithms implemented compute    *
 * the latter definitions. We have taken this into account in our functions    *
 * and return the values corresponding to the equations above.                 *
 *******************************************************************************
 *                              DEFINED FUNCTIONS                              *
 *******************************************************************************
 *  Fresnel_Sine_Taylor_to_Asymptotic / Fresnel_Cosine_Taylor_to_Asymptotic:   *
 *      This uses the standard Taylor expansion for small inputs (|x|<=3), and *
 *      asymptotic expansions for large input (|x|>3). The Taylor Series are:  *
 *                                                                             *
 *                  ___                                                      *
 *                  \      (-1)^n x^(4n+3)/                                    *
 *        S(x)=     /__                /  (4n+3)(2n+1)!                      *
 *                  n = 0                                                      *
 *                                                                             *
 *      This can be obtained by substituting the Taylor expansions for         *
 *      sin(x^2) and cos(x^2), respectively, and integrating term by term.     *
 *      The interval [0,3] is broken into three pieces, and the appropriate    *
 *      number of terms in the expansion are used accordingly to obtain at     *
 *      least 8 decimals of accuracy.                                          *
 *                                                                             *
 *      The asymptotic expansions can be obtained by iteratively using         *
 *      integration by parts. The Asymptotic expansions diverge for all x, but *
 *      by cutting off the expansion at a particular N, one finds a very good  *
 *      approximations. The series are given by:                               *
 *                                                                             *
 *          a_n(x) = (4n+2)! / (2^(4n+3) (2n+1)! x^(4n+3))                     *
 *          b_n(x) = (4n)! / (2^(4n+1) (2n)! x^(4n+1))                         *
 *                                                                             *
 *                             ___                                           *
 *                             \                                               *
 *        S(x) = sqrt(pi/i) -  /__  (-1)^n (a_n(x)sin(x^2) + b_n(x)cos(x^2)) *
 *                             n = 0                                           *
 *                                                                             *
 *      The error in the asympotic series goes like |a_N(x)|+|b_N(x)|.         *
 *      For large x, and appropriate N, this can be made incredibly small.     *
 *******************************************************************************
 *  Fresnel_Sine_While_to_Asymptotic / Fresnel_Cosine_While_to_Asymptotic:     *
 *      Similar to Fresnel_Sine_Taylor_to_Asymptotic, but a while loop is used *
 *      during the Taylor approximation, stopping once the error is 1.0e-8.    *
 *******************************************************************************
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Three,                               *
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Four,                                *
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Six,                                 *
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Eight,                               *
 *      Computes fresnel integrals to specified precision using a rational     *
 *      approximation as computed by Mark A. Heald.                            *
 *      See Rational Approximations for the Fresnel Integrals, Mark A. Heald,  *
 *      Mathematics of Computation, Vol. 44, No. 170 (Apr., 1985), pp. 459-461 *
 *******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                                *
 *  Date:       Febuary 26, 2019                                               *
 ******************************************************************************/
#include <math.h>

/*  Various coefficients and constants defined here.                          */
#include "_math_constants.h"

/*  Various Fresnel integral functions declared here.                         */
#include "special_functions.h"

/*******************************************************************************
 *------------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.   *
 * They are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.                *
 ******************************************************************************/

float Fresnel_Sin_Float(float x)
{
    /* Variables for S(x) and powers of x, respectively. */
    float sx;
    float arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 9.0)
    {
        x *= arg;
        arg *= arg;
        sx = arg * FRESNEL_SINE_TAYLOR_16 + FRESNEL_SINE_TAYLOR_15;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_14;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_13;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_12;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_11;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_10;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_09;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_08;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_07;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_06;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_05;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_04;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_03;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_02;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_01;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_00;
        return sx*x;
    }
    else if (arg < 1.0e16) {
        float sinarg, cosarg, cos_x_squared, sin_x_squared;
        cos_x_squared = cosf(arg);
        sin_x_squared = sinf(arg);

        arg = 1.0/arg;
        cos_x_squared *= arg;
        arg *= arg;
        sin_x_squared *= arg;

        cosarg  = arg * FRESNEL_SINE_ASYM_08 + FRESNEL_SINE_ASYM_06;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_04;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_02;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_00;
        cosarg *= cos_x_squared;

        sinarg  = arg * FRESNEL_SINE_ASYM_09 + FRESNEL_SINE_ASYM_07;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_05;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_03;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_01;
        sinarg *= sin_x_squared;

        sx = cosarg + sinarg;
        sx *= x;

        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return sx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    /*  For large values, return the limit of S(x) as x -> +/- infinity.      */
    else return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
}

double Fresnel_Sin_Double(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double sx;
    double arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 11.68){
        x *= arg;
        arg *= arg;
        sx = arg * FRESNEL_SINE_TAYLOR_22 + FRESNEL_SINE_TAYLOR_21;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_20;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_19;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_18;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_17;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_16;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_15;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_14;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_13;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_12;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_11;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_10;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_09;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_08;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_07;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_06;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_05;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_04;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_03;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_02;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_01;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_00;
        return sx*x;
    }
    else if (arg < 1.0e16) {
        double sinarg, cosarg, cos_x_squared, sin_x_squared;
        cos_x_squared = cos(arg);
        sin_x_squared = sin(arg);

        arg = 1.0/arg;
        cos_x_squared *= arg;
        arg *= arg;
        sin_x_squared *= arg;

        cosarg  = arg * FRESNEL_SINE_ASYM_08 + FRESNEL_SINE_ASYM_06;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_04;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_02;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_00;
        cosarg *= cos_x_squared;

        sinarg  = arg * FRESNEL_SINE_ASYM_09 + FRESNEL_SINE_ASYM_07;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_05;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_03;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_01;
        sinarg *= sin_x_squared;

        sx = cosarg + sinarg;
        sx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return sx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    /*  For large values, return the limit of S(x) as x -> +/- infinity.      */
    else return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
}

long double Fresnel_Sin_Long_Double(long double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    long double sx;
    long double arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 11.68)
    {
        x *= arg;
        arg *= arg;
        sx = arg * FRESNEL_SINE_TAYLOR_24 + FRESNEL_SINE_TAYLOR_23;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_22;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_21;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_20;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_19;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_18;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_17;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_16;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_15;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_14;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_13;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_12;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_11;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_10;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_09;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_08;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_07;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_06;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_05;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_04;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_03;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_02;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_01;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_00;
        return sx*x;
    }
    else if (arg < 1.0e16) {
        long double sinarg, cosarg, cos_x_squared, sin_x_squared;
        cos_x_squared = cosl(arg);
        sin_x_squared = sinl(arg);

        arg = 1.0/arg;
        cos_x_squared *= arg;
        arg *= arg;
        sin_x_squared *= arg;

        cosarg  = arg * FRESNEL_SINE_ASYM_08 + FRESNEL_SINE_ASYM_06;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_04;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_02;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_00;
        cosarg *= cos_x_squared;

        sinarg  = arg * FRESNEL_SINE_ASYM_09 + FRESNEL_SINE_ASYM_07;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_05;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_03;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_01;
        sinarg *= sin_x_squared;

        sx = cosarg + sinarg;
        sx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return sx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    /* For large values, return the limit of S(x) as x -> +/- infinity.       */
    else return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
}

/*  For all integer types, convert to double and compute.                     */
RSSRINGOCCSNonFloatInputForFloatOutput(Fresnel_Sin);
