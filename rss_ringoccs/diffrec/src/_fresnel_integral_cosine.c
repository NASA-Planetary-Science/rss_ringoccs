/*******************************************************************************
 *                             Fresnel Cosine                                  *
 *******************************************************************************
 *  This C program contains several algorithms for the computation of the      *
 *  Fresnel Cosine integral. This is meant to test the accuracy and efficiency *
 *  of the various algorithms. This file is not included in the setup.py file, *
 *  and thus will not be compiled. The fresnel sine and cosine functions that  *
 *  are used in the _special_functions.c file are the functions that had the   *
 *  best combination of speed and accuracy. This file is kept so that users    *
 *  may experiment with the various known algorithms.                          *
 *******************************************************************************
 *  We define the Fresnel Cosine Integrals as follows:                         *
 *                 x                                                           *
 *                 -                                                           *
 *                | |                                                          *
 *       C(x) =   |   cos(t^2) dt                                              *
 *              | |                                                            *
 *               -                                                             *
 *               0                                                             *
 *******************************************************************************
 *  It is very common for a pi/2 to be placed inside the cosine term,          *
 *  and thus in translating one would need to scale x by sqrt(2/pi) and scale  *
 *  the results by sqrt(pi/2). Several of the algorithms implemented compute   *
 *  the latter definition. We have taken this into account in our functions    *
 *  and return the values corresponding to the equations above.                *
 *******************************************************************************
 *                              DEFINED FUNCTIONS                              *
 *******************************************************************************
 *  Fresnel_Sine_Taylor_to_Asymptotic:                                         *
 *      This uses the standard Taylor expansion for small inputs (|x|<=4), and *
 *      asymptotic expansions for large input (|x|>4). The Taylor Series are:  *
 *                  ___                                                      *
 *                  \      (-1)^n x^(4n+1)/                                    *
 *        C(x)=     /__                /  (4n+1)(2n)!                        *
 *                  n = 0                                                      *
 *                                                                             *
 *      This can be obtained by substituting the Taylor expansions for         *
 *      sin(x^2) and cos(x^2), respectively, and integrating term by term.     *
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
 *        C(x) = sqrt(pi/i) +  /__  (-1)^n (b_n(x)sin(x^2) - a_n(x)cos(x^2)) *
 *                             n = 0                                           *
 *                                                                             *
 *      The error in the asympotic series goes like |a_N(x)|+|b_N(x)|.         *
 *      For large x, and appropriate N, this can be made incredibly small.     *
 *      For |x|>10, and choosing N=6, the max error is 1.e-12. For values      *
 *      near |x|=4 the error is 1.e-6.                                         *
 *******************************************************************************
 *  Fresnel_Sine_While_to_Asymptotic:                                          *
 *      Similar to Fresnel_Sine_Taylor_to_Asymptotic, but a while loop is used *
 *      during the Taylor approximation, stopping once the error is 1.0e-8.    *
 *******************************************************************************
 *  Fresnel_Cosine_Heald_Rational_EPS_Minus_Three,                             *
 *  Fresnel_Cosine_Heald_Rational_EPS_Minus_Three,                             *
 *  Fresnel_Cosine_Heald_Rational_EPS_Minus_Six,                               *
 *  Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight:                             *
 *      Computes fresnel integrals to specified precision using a rational     *
 *      approximation as computed by Mark A. Heald.                            *
 *      See Rational Approximations for the Fresnel Integrals, Mark A. Heald,  *
 *      Mathematics of Computation, Vol. 44, No. 170 (Apr., 1985), pp. 459-461 *
 *******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                                *
 *  Date:       August 30, 2019                                                *
 ******************************************************************************/
/*  Various coefficients and constants defined here.                          */
#include "_math_constants.h"

/*  Various Fresnel integral functions declared here.                         */
#include "special_functions.h"

float Fresnel_Cos_Float(float x)
{
    /* Variables for S(x) and powers of x, respectively. */
    float cx, arg;
    arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        arg *= arg;
        cx = arg * FRESNEL_COSINE_TAYLOR_15 + FRESNEL_COSINE_TAYLOR_14;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_13;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_12;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_11;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_10;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_09;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_08;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_07;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_06;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_05;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_04;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_03;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_02;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_01;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_00;
        return cx*x;
    }
    else if (arg < 1.0e16) {
        float sinarg, cosarg, cos_x_squared, sin_x_squared;
        cos_x_squared = cosf(arg);
        sin_x_squared = sinf(arg);

        arg = 1.0/arg;
        sin_x_squared *= arg;
        arg *= arg;
        cos_x_squared *= arg;

        sinarg  = arg * FRESNEL_COSINE_ASYM_08 + FRESNEL_COSINE_ASYM_06;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_04;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_02;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_00;
        sinarg *= sin_x_squared;

        cosarg  = arg * FRESNEL_COSINE_ASYM_09 + FRESNEL_COSINE_ASYM_07;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_05;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_03;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_01;
        cosarg *= cos_x_squared;

        cx = cosarg + sinarg;
        cx *= x;

        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return cx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    else {
        /* For large values, return the limit of S(x) as x -> +/- infinity. */
        return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
}

double Fresnel_Cos_Double(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double cx, arg;
    arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 13.19){
        arg *= arg;
        cx = arg * FRESNEL_COSINE_TAYLOR_22 + FRESNEL_COSINE_TAYLOR_21;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_20;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_19;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_18;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_17;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_16;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_15;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_14;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_13;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_12;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_11;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_10;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_09;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_08;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_07;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_06;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_05;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_04;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_03;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_02;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_01;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_00;
        return cx*x;
    }
    else if (arg < 1.0e16) {
        double sinarg, cosarg, cos_x_squared, sin_x_squared;
        cos_x_squared = cos(arg);
        sin_x_squared = sin(arg);

        arg = 1.0/arg;
        sin_x_squared *= arg;
        arg *= arg;
        cos_x_squared *= arg;

        sinarg  = arg * FRESNEL_COSINE_ASYM_08 + FRESNEL_COSINE_ASYM_06;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_04;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_02;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_00;
        sinarg *= sin_x_squared;

        cosarg  = arg * FRESNEL_COSINE_ASYM_09 + FRESNEL_COSINE_ASYM_07;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_05;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_03;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_01;
        cosarg *= cos_x_squared;

        cx = cosarg + sinarg;
        cx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return cx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    else {
        /* For large values, return the limit of S(x) as x -> +/- infinity. */
        return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
}

long double Fresnel_Cos_Long_Double(long double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    long double cx, arg;
    arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 16.24){
        arg *= arg;
        cx = arg * FRESNEL_COSINE_TAYLOR_26 + FRESNEL_COSINE_TAYLOR_25;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_24;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_23;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_22;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_21;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_20;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_19;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_18;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_17;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_16;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_15;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_14;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_13;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_12;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_11;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_10;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_09;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_08;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_07;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_06;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_05;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_04;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_03;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_02;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_01;
        cx = arg * cx + FRESNEL_COSINE_TAYLOR_00;
        return cx*x;
    }
    else if (arg < 1.0e16) {
        long double sinarg, cosarg, cos_x_squared, sin_x_squared;
        cos_x_squared = cosl(arg);
        sin_x_squared = sinl(arg);

        arg = 1.0/arg;
        sin_x_squared *= arg;
        arg *= arg;
        cos_x_squared *= arg;

        sinarg  = arg * FRESNEL_COSINE_ASYM_08 + FRESNEL_COSINE_ASYM_06;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_04;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_02;
        sinarg  = arg * sinarg + FRESNEL_COSINE_ASYM_00;
        sinarg *= sin_x_squared;

        cosarg  = arg * FRESNEL_COSINE_ASYM_09 + FRESNEL_COSINE_ASYM_07;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_05;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_03;
        cosarg  = arg * cosarg + FRESNEL_COSINE_ASYM_01;
        cosarg *= cos_x_squared;

        cx = cosarg + sinarg;
        cx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return cx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    else {
        /* For large values, return the limit of S(x) as x -> +/- infinity. */
        return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
}

/*  For all integer types, convert to double and compute.                     */
RSSRINGOCCSNonFloatInputForFloatOutput(Fresnel_Cos);
