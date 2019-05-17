/*******************************************************************************
 *                             Fresnel Cosine                                  *
 *******************************************************************************
 * This C program contains several algorithms for the computation of the       *
 * Fresnel Cosine integral. This is meant to test the accuracy and efficiency  *
 * of the various algorithms. This file is not included in the setup.py file,  *
 * and thus will not be compiled. The fresnel sine and cosine functions that   *
 * are used in the _special_functions.c file are the functions that had the    *
 * best combination of speed and accuracy. This file is kept so that users may *
 * experiment with the various known algorithms.                               *
 *******************************************************************************
 * We define the Fresnel Cosine Integrals as follows:                          *
 *           x                                                                 *
 *           -                                                                 *
 *          | |                                                                *
 * C(x) =   |   cos(t^2) dt                                                    *
 *        | |                                                                  *
 *         -                                                                   *
 *         0                                                                   *
 *******************************************************************************
 * It is very common for a pi/2 to be placed inside the cosine term,           *
 * and thus in translating one would need to scale x by sqrt(2/pi) and scale   *
 * the results by sqrt(pi/2). Several of the algorithms implemented compute    *
 * the latter definition. We have taken this into account in our functions     *
 * and return the values corresponding to the equations above.                 *
 *******************************************************************************
 *                              DEFINED FUNCTIONS                              *
 *******************************************************************************
 *  Fresnel_Sine_Taylor_to_Asymptotic:                                         *
 *      This uses the standard Taylor expansion for small inputs (|x|<=4), and *
 *      asymptotic expansions for large input (|x|>4). The Taylor Series are:  *
 *                  _____                                                      *
 *                  \      (-1)^n x^(4n+1)/                                    *
 *        C(x)=     /____                /  (4n+1)(2n)!                        *
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
 *                             _____                                           *
 *                             \                                               *
 *        C(x) = sqrt(pi/i) +  /____  (-1)^n (b_n(x)sin(x^2) - a_n(x)cos(x^2)) *
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
 *  Date:       Febuary 26, 2019                                               *
 ******************************************************************************/
#ifndef RSS_RINGOCCS_FRESNEL_COSINE_H
#define RSS_RINGOCCS_FRESNEL_COSINE_H

#include <math.h>
#include <complex.h>

/*  Various coefficients and constants defined here.    */
#include "__math_constants.h"

/*******************************************************************************
 *------------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.   *
 * They are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.                *
 ******************************************************************************/

double Fresnel_Cosine_Taylor_to_Asymptotic_Func(double x)
{
    /***************************************************************************
     * This is the primary function for compute the Fresnel_Cosine integral.   *
     * It offers the best combination of speed and accurary. For small         *
     * arguments a Taylor Series is computed to 8 decimals of accuracy, 35 or  *
     * more for arguments less than one. In this region accuracy is limited by *
     * double precision round-off error, and not by the Taylor series.         *
     *                                                                         *
     * For Large values the asymptotic expansion is used. This error for this  *
     * is e-6 in the vicinity of |x| = 4, but quickly decays to zero for       *
     * larger values. For |x|>10, double precision accuracy is achieved.       *
     *                                                                         *
     * For extremely large values, return the limit of the Fresnel Cosine      *
     * function. This is +/- SQRT(pi/8). The error is within double precision  *
     * for when the argument is greater than 1.e16.                            *
     **************************************************************************/

    /* Variables for S(x) and powers of x, respectively. */
    double cx, arg;
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
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        sinarg *= arg;
        arg *= arg;
        cosarg *= arg;

        sinarg *= FRESNEL_COSINE_ASYM_00 + arg*(
                    FRESNEL_COSINE_ASYM_02 + arg*(
                        FRESNEL_COSINE_ASYM_04 + FRESNEL_COSINE_ASYM_06*arg
                    )
                );
        cosarg *= FRESNEL_SINE_ASYM_01 + arg*(
                    FRESNEL_SINE_ASYM_03 + arg*(
                        FRESNEL_SINE_ASYM_05 + FRESNEL_SINE_ASYM_07*arg
                    )
                );

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

double Fresnel_Cosine_While_to_Asymptotic_Func(double x)
{
    double FRESNEL_COSINE_TAYLOR_COEFFICIENTS[27] = {
        FRESNEL_COSINE_TAYLOR_00, FRESNEL_COSINE_TAYLOR_01,
        FRESNEL_COSINE_TAYLOR_02, FRESNEL_COSINE_TAYLOR_03,
        FRESNEL_COSINE_TAYLOR_04, FRESNEL_COSINE_TAYLOR_05,
        FRESNEL_COSINE_TAYLOR_06, FRESNEL_COSINE_TAYLOR_07,
        FRESNEL_COSINE_TAYLOR_08, FRESNEL_COSINE_TAYLOR_09,
        FRESNEL_COSINE_TAYLOR_10, FRESNEL_COSINE_TAYLOR_11,
        FRESNEL_COSINE_TAYLOR_12, FRESNEL_COSINE_TAYLOR_13,
        FRESNEL_COSINE_TAYLOR_14, FRESNEL_COSINE_TAYLOR_15,
        FRESNEL_COSINE_TAYLOR_16, FRESNEL_COSINE_TAYLOR_17,
        FRESNEL_COSINE_TAYLOR_18, FRESNEL_COSINE_TAYLOR_19,
        FRESNEL_COSINE_TAYLOR_20, FRESNEL_COSINE_TAYLOR_21,
        FRESNEL_COSINE_TAYLOR_22, FRESNEL_COSINE_TAYLOR_23,
        FRESNEL_COSINE_TAYLOR_24, FRESNEL_COSINE_TAYLOR_25,
        FRESNEL_COSINE_TAYLOR_26
    };

    /* Variables for S(x) and powers of x, respectively. */
    double cx;
    double arg = x*x;
    double x4 = arg*arg;
    double EPS = 1.0e-8;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        int i = 0;
        double term = arg*FRESNEL_COSINE_TAYLOR_COEFFICIENTS[0];
        cx = term;
        while (term > EPS){
            /* Odd terms have a negative coefficients.
               Compute two terms at a time to compare with EPS. */
            i += 1;
            arg *= x4;
            term = arg*FRESNEL_COSINE_TAYLOR_COEFFICIENTS[i];
            cx += term;

            // Compute even term.
            i += 1;
            arg *= x4;
            term = arg*FRESNEL_COSINE_TAYLOR_COEFFICIENTS[i];
            cx += term;
        }
        return x*cx;
    }
    else if (arg < 1.0e16) {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        sinarg *= arg;
        arg *= arg;
        cosarg *= arg;

        sinarg *= FRESNEL_COSINE_ASYM_00 + arg*(
                    FRESNEL_COSINE_ASYM_02 + arg*(
                        FRESNEL_COSINE_ASYM_04 + FRESNEL_COSINE_ASYM_06*arg
                    )
                );
        cosarg *= FRESNEL_SINE_ASYM_01 + arg*(
                    FRESNEL_SINE_ASYM_03 + arg*(
                        FRESNEL_SINE_ASYM_05 + FRESNEL_SINE_ASYM_07*arg
                    )
                );

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

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Three(double x)
{
    double A, R, a, b, c, d, sgn_x;
    sgn_x = (x>0)-(x<0);
    x *= SQRT_2_BY_PI*sgn_x;

    /* Compute the Numerator of the A_jk Function.      */
    a = FRESNEL_HEALD_RATIONAL_EPS_3_A00;

    /* Compute the Denominator of the A_jk Function.    */
    b = FRESNEL_HEALD_RATIONAL_EPS_3_B03*x + FRESNEL_HEALD_RATIONAL_EPS_3_B02;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_3_B01;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_3_B00;

    /* Compute the Numerator of the R_lm Function.      */
    c = FRESNEL_HEALD_RATIONAL_EPS_3_C01*x + FRESNEL_HEALD_RATIONAL_EPS_3_C00;

    /* Compute the Denominator of the R_lm Function.    */
    d = FRESNEL_HEALD_RATIONAL_EPS_3_D02*x + FRESNEL_HEALD_RATIONAL_EPS_3_D01;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_3_D00;

    A = a/b-x*x;
    A *= PI_BY_TWO;
    R = c/d;
    R *= SQRT_PI_BY_2;

    return sgn_x*(SQRT_PI_BY_8 - R*sin(A));
}

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Four(double x)
{
    double A, R, a, b, c, d, sgn_x;
    sgn_x = (x>0)-(x<0);
    x *= SQRT_2_BY_PI*sgn_x;

    /* Compute the Numerator of the A_jk Function.      */
    a = FRESNEL_HEALD_RATIONAL_EPS_4_A01*x + FRESNEL_HEALD_RATIONAL_EPS_4_A00;

    /* Compute the Denominator of the A_jk Function.    */
    b = FRESNEL_HEALD_RATIONAL_EPS_4_B03*x + FRESNEL_HEALD_RATIONAL_EPS_4_B02;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_4_B01;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_4_B00;

    /* Compute the Numerator of the R_lm Function.      */
    c = FRESNEL_HEALD_RATIONAL_EPS_4_C02*x + FRESNEL_HEALD_RATIONAL_EPS_4_C01;
    c = x*c + FRESNEL_HEALD_RATIONAL_EPS_4_C00;

    /* Compute the Denominator of the R_lm Function.    */
    d = FRESNEL_HEALD_RATIONAL_EPS_4_D03*x + FRESNEL_HEALD_RATIONAL_EPS_4_D02;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_4_D01;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_4_D00;

    A = a/b-x*x;
    A *= PI_BY_TWO;
    R = c/d;
    R *= SQRT_PI_BY_2;

    return sgn_x*(SQRT_PI_BY_8 - R*sin(A));
}

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Six(double x)
{
    double A, R, a, b, c, d, sgn_x;
    sgn_x = (x>0)-(x<0);
    x *= SQRT_2_BY_PI*sgn_x;

    /* Compute the Numerator of the A_jk Function.      */
    a = FRESNEL_HEALD_RATIONAL_EPS_6_A02*x + FRESNEL_HEALD_RATIONAL_EPS_6_A01;
    a = a*x + FRESNEL_HEALD_RATIONAL_EPS_6_A00;

    /* Compute the Denominator of the A_jk Function.    */
    b = FRESNEL_HEALD_RATIONAL_EPS_6_B04*x + FRESNEL_HEALD_RATIONAL_EPS_6_B03;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_6_B02;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_6_B01;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_6_B00;

    /* Compute the Numerator of the R_lm Function.      */
    c = FRESNEL_HEALD_RATIONAL_EPS_6_C03*x + FRESNEL_HEALD_RATIONAL_EPS_6_C02;
    c = x*c + FRESNEL_HEALD_RATIONAL_EPS_6_C01;
    c = x*c + FRESNEL_HEALD_RATIONAL_EPS_6_C00;

    /* Compute the Denominator of the R_lm Function.    */
    d = FRESNEL_HEALD_RATIONAL_EPS_6_D04*x + FRESNEL_HEALD_RATIONAL_EPS_6_D03;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_6_D02;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_6_D01;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_6_D00;

    A = a/b-x*x;
    A *= PI_BY_TWO;
    R = c/d;
    R *= SQRT_PI_BY_2;

    return sgn_x*(SQRT_PI_BY_8 - R*sin(A));
}

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight(double x)
{
    double A, R, a, b, c, d, sgn_x;
    sgn_x = (x>0)-(x<0);
    x *= SQRT_2_BY_PI*sgn_x;

    /* Compute the Numerator of the A_jk Function.      */
    a = FRESNEL_HEALD_RATIONAL_EPS_8_A04*x + FRESNEL_HEALD_RATIONAL_EPS_8_A03;
    a = a*x + FRESNEL_HEALD_RATIONAL_EPS_8_A02;
    a = a*x + FRESNEL_HEALD_RATIONAL_EPS_8_A01;
    a = a*x + FRESNEL_HEALD_RATIONAL_EPS_8_A00;

    /* Compute the Denominator of the A_jk Function.    */
    b = FRESNEL_HEALD_RATIONAL_EPS_8_B06*x + FRESNEL_HEALD_RATIONAL_EPS_8_B05;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_8_B04;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_8_B03;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_8_B02;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_8_B01;
    b = x*b + FRESNEL_HEALD_RATIONAL_EPS_8_B00;

    /* Compute the Numerator of the R_lm Function.      */
    c = FRESNEL_HEALD_RATIONAL_EPS_8_C05*x + FRESNEL_HEALD_RATIONAL_EPS_8_C04;
    c = x*c + FRESNEL_HEALD_RATIONAL_EPS_8_C03;
    c = x*c + FRESNEL_HEALD_RATIONAL_EPS_8_C02;
    c = x*c + FRESNEL_HEALD_RATIONAL_EPS_8_C01;
    c = x*c + FRESNEL_HEALD_RATIONAL_EPS_8_C00;

    /* Compute the Denominator of the R_lm Function.    */
    d = FRESNEL_HEALD_RATIONAL_EPS_8_D06*x + FRESNEL_HEALD_RATIONAL_EPS_8_D05;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_8_D04;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_8_D03;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_8_D02;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_8_D01;
    d = x*d + FRESNEL_HEALD_RATIONAL_EPS_8_D00;

    A = a/b-x*x;
    A *= PI_BY_TWO;
    R = c/d;
    R *= SQRT_PI_BY_2;

    return sgn_x*(SQRT_PI_BY_8 - R*sin(A));
}

#endif