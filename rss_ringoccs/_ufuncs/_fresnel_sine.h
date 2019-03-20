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
 *                  _____                                                      *
 *                  \      (-1)^n x^(4n+3)/                                    *
 *        S(x)=     /____                /  (4n+3)(2n+1)!                      *
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
 *                             _____                                           *
 *                             \                                               *
 *        S(x) = sqrt(pi/i) -  /____  (-1)^n (a_n(x)sin(x^2) + b_n(x)cos(x^2)) *
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
#include <complex.h>

/* Define Miscellaneous Constants. */
#define SQRT_PI_BY_8 0.6266570686577501
#define SQRT_PI_BY_2 1.2533141373155001
#define SQRT_2_BY_PI 0.7978845608028654
#define PI_BY_TWO 1.5707963267948966
#define SQRT_2 1.4142135623730951

/* Define Coefficients for the Fresnel Sine Taylor Expansion. */
#define FRESNEL_SINE_TAYLOR_00 0.3333333333333333
#define FRESNEL_SINE_TAYLOR_01 -0.023809523809523808
#define FRESNEL_SINE_TAYLOR_02 0.0007575757575757576
#define FRESNEL_SINE_TAYLOR_03 -1.3227513227513228e-05
#define FRESNEL_SINE_TAYLOR_04 1.4503852223150468e-07
#define FRESNEL_SINE_TAYLOR_05 -1.0892221037148573e-09
#define FRESNEL_SINE_TAYLOR_06 5.9477940136376354e-12
#define FRESNEL_SINE_TAYLOR_07 -2.466827010264457e-14
#define FRESNEL_SINE_TAYLOR_08 8.032735012415773e-17
#define FRESNEL_SINE_TAYLOR_09 -2.107855191442136e-19
#define FRESNEL_SINE_TAYLOR_10 4.5518467589282e-22
#define FRESNEL_SINE_TAYLOR_11 -8.230149299214221e-25
#define FRESNEL_SINE_TAYLOR_12 1.2641078988989164e-27
#define FRESNEL_SINE_TAYLOR_13 -1.669761793417372e-30
#define FRESNEL_SINE_TAYLOR_14 1.9169428621097826e-33
#define FRESNEL_SINE_TAYLOR_15 -1.9303572088151077e-36
#define FRESNEL_SINE_TAYLOR_16 1.7188560628017835e-39
#define FRESNEL_SINE_TAYLOR_17 -1.3630412617791397e-42
#define FRESNEL_SINE_TAYLOR_18 9.687280238870763e-46
#define FRESNEL_SINE_TAYLOR_19 -6.205657919637396e-49
#define FRESNEL_SINE_TAYLOR_20 3.601579309810126e-52
#define FRESNEL_SINE_TAYLOR_21 -1.9025412272898796e-55
#define FRESNEL_SINE_TAYLOR_22 9.186429502398686e-59
#define FRESNEL_SINE_TAYLOR_23 -4.070135277853256e-62
#define FRESNEL_SINE_TAYLOR_24 1.66058051345109e-65
#define FRESNEL_SINE_TAYLOR_25 6.259184116948712e-69
#define FRESNEL_SINE_TAYLOR_26 2.1862104229538858e-72

/* Define Coefficients for the Fresnel Sine Asymptotic Expansion. */
#define FRESNEL_SINE_ASYM_00 -0.5
#define FRESNEL_SINE_ASYM_01 -0.25
#define FRESNEL_SINE_ASYM_02 0.375
#define FRESNEL_SINE_ASYM_03 0.9375
#define FRESNEL_SINE_ASYM_04 -3.281250
#define FRESNEL_SINE_ASYM_05 -14.765625
#define FRESNEL_SINE_ASYM_06 81.210938
#define FRESNEL_SINE_ASYM_07 527.87109375

/*******************************************************************************
 * Define Coefficients Used for the Rational Approximation of the              *
 * Fresnel Integrals using approximations from Mark. A. Heald.                 *
 * See Rational Approximations for the Fresnel Integrals,                      *
 * Mathematics of Computation, Vol. 44, No. 170 (Apr., 1985), pp. 459-461      *
 *                                                                             *
 * Heald defines the Fresnel Integrals as the integral of                      *
 * sin(pi/2 x^2) and cos(pi/2 x^2), whereas we adopt the definition of         *
 * the integral of sin(x^2) and cos(x^2). As such, a scale factor of           *
 * sqrt(2/pi) is multiplied to our coefficients, and a scale factor of         *
 * sqrt(pi/2) is multiplied to the final output.                               *
 ******************************************************************************/

/* Coefficients for up to 3 significant digits. */
#define FRESNEL_HEALD_RATIONAL_EPS_3_A00 1.00000

#define FRESNEL_HEALD_RATIONAL_EPS_3_B00 2.000000
#define FRESNEL_HEALD_RATIONAL_EPS_3_B01 2.524
#define FRESNEL_HEALD_RATIONAL_EPS_3_B02 1.886
#define FRESNEL_HEALD_RATIONAL_EPS_3_B03 0.803

#define FRESNEL_HEALD_RATIONAL_EPS_3_C00 1.00000
#define FRESNEL_HEALD_RATIONAL_EPS_3_C01 0.506

#define FRESNEL_HEALD_RATIONAL_EPS_3_D00 1.4142135623730951
#define FRESNEL_HEALD_RATIONAL_EPS_3_D01 2.054
#define FRESNEL_HEALD_RATIONAL_EPS_3_D02 1.79

/* Coefficients for up to 4 significant digits. */
#define FRESNEL_HEALD_RATIONAL_EPS_4_A00 1.00000
#define FRESNEL_HEALD_RATIONAL_EPS_4_A01 0.1765

#define FRESNEL_HEALD_RATIONAL_EPS_4_B00 2.0000
#define FRESNEL_HEALD_RATIONAL_EPS_4_B01 2.915
#define FRESNEL_HEALD_RATIONAL_EPS_4_B02 2.079
#define FRESNEL_HEALD_RATIONAL_EPS_4_B03 1.519

#define FRESNEL_HEALD_RATIONAL_EPS_4_C00 1.00000
#define FRESNEL_HEALD_RATIONAL_EPS_4_C01 0.5083
#define FRESNEL_HEALD_RATIONAL_EPS_4_C02 0.3569

#define FRESNEL_HEALD_RATIONAL_EPS_4_D00 1.4142135623730951
#define FRESNEL_HEALD_RATIONAL_EPS_4_D01 2.1416
#define FRESNEL_HEALD_RATIONAL_EPS_4_D02 1.8515
#define FRESNEL_HEALD_RATIONAL_EPS_4_D03 1.1021

/* Coefficients for up to 6 significant digits. */
#define FRESNEL_HEALD_RATIONAL_EPS_6_A00 1.00000
#define FRESNEL_HEALD_RATIONAL_EPS_6_A01 0.08218
#define FRESNEL_HEALD_RATIONAL_EPS_6_A02 0.15108

#define FRESNEL_HEALD_RATIONAL_EPS_6_B00 2.0000
#define FRESNEL_HEALD_RATIONAL_EPS_6_B01 2.7097
#define FRESNEL_HEALD_RATIONAL_EPS_6_B02 2.3185
#define FRESNEL_HEALD_RATIONAL_EPS_6_B03 1.2389
#define FRESNEL_HEALD_RATIONAL_EPS_6_B04 0.6561

#define FRESNEL_HEALD_RATIONAL_EPS_6_C00 1.00000
#define FRESNEL_HEALD_RATIONAL_EPS_6_C01 0.60427
#define FRESNEL_HEALD_RATIONAL_EPS_6_C02 0.41159
#define FRESNEL_HEALD_RATIONAL_EPS_6_C03 0.19170

#define FRESNEL_HEALD_RATIONAL_EPS_6_D00 1.4142135623730951
#define FRESNEL_HEALD_RATIONAL_EPS_6_D01 2.26794
#define FRESNEL_HEALD_RATIONAL_EPS_6_D02 2.15594
#define FRESNEL_HEALD_RATIONAL_EPS_6_D03 1.26057
#define FRESNEL_HEALD_RATIONAL_EPS_6_D04 0.60353

/* Coefficients for up to 8 significant digits. */
#define FRESNEL_HEALD_RATIONAL_EPS_8_A00 1.0000000
#define FRESNEL_HEALD_RATIONAL_EPS_8_A01 0.1945161
#define FRESNEL_HEALD_RATIONAL_EPS_8_A02 0.2363641
#define FRESNEL_HEALD_RATIONAL_EPS_8_A03 0.0683240
#define FRESNEL_HEALD_RATIONAL_EPS_8_A04 0.0241212

#define FRESNEL_HEALD_RATIONAL_EPS_8_B00 2.0000000
#define FRESNEL_HEALD_RATIONAL_EPS_8_B01 2.9355041
#define FRESNEL_HEALD_RATIONAL_EPS_8_B02 2.7570460
#define FRESNEL_HEALD_RATIONAL_EPS_8_B03 1.8757210
#define FRESNEL_HEALD_RATIONAL_EPS_8_B04 0.9781130
#define FRESNEL_HEALD_RATIONAL_EPS_8_B05 0.3566810
#define FRESNEL_HEALD_RATIONAL_EPS_8_B06 0.1182470

#define FRESNEL_HEALD_RATIONAL_EPS_8_C00 1.0000000
#define FRESNEL_HEALD_RATIONAL_EPS_8_C01 0.7769507
#define FRESNEL_HEALD_RATIONAL_EPS_8_C02 0.6460117
#define FRESNEL_HEALD_RATIONAL_EPS_8_C03 0.3460509
#define FRESNEL_HEALD_RATIONAL_EPS_8_C04 0.1339259
#define FRESNEL_HEALD_RATIONAL_EPS_8_C05 0.0433995

#define FRESNEL_HEALD_RATIONAL_EPS_8_D00 1.4142135623730951
#define FRESNEL_HEALD_RATIONAL_EPS_8_D01 2.5129806
#define FRESNEL_HEALD_RATIONAL_EPS_8_D02 2.7196741
#define FRESNEL_HEALD_RATIONAL_EPS_8_D03 1.9840524
#define FRESNEL_HEALD_RATIONAL_EPS_8_D04 1.0917325
#define FRESNEL_HEALD_RATIONAL_EPS_8_D05 0.4205217
#define FRESNEL_HEALD_RATIONAL_EPS_8_D06 0.13634704

/*******************************************************************************
 *------------------------------DEFINE C FUNCTIONS-----------------------------*
 *******************************************************************************
 * These are functions written in pure C without the use of the Numpy-C API.   *
 * The are used to define various special functions. They will be wrapped in   *
 * a form that is useable with the Python interpreter later on.                *
 ******************************************************************************/
double Fresnel_Sine_Taylor_to_Asymptotic_Func(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double sx;
    double arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 9.0){
        x *= arg;
        arg *= arg;
        sx = arg * FRESNEL_SINE_TAYLOR_15 + FRESNEL_SINE_TAYLOR_14;
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
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        cosarg *= arg;
        arg *= arg;
        sinarg *= arg;

        cosarg *= FRESNEL_SINE_ASYM_00 + arg*(
                    FRESNEL_SINE_ASYM_02 + arg*(
                        FRESNEL_SINE_ASYM_04 + FRESNEL_SINE_ASYM_06*arg
                    )
                );
        sinarg *= FRESNEL_SINE_ASYM_01 + arg*(
                    FRESNEL_SINE_ASYM_03 + arg*(
                        FRESNEL_SINE_ASYM_05 + FRESNEL_SINE_ASYM_07*arg
                    )
                );

        sx = cosarg + sinarg;
        sx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return sx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    else {
        /* For large values, return the limit of S(x) as x -> +/- infinity. */
        return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
}

double Fresnel_Sine_While_to_Asymptotic_Func(double x)
{
    double FRESNEL_SINE_TAYLOR_COEFFICIENTS[27] = {
        FRESNEL_SINE_TAYLOR_00, FRESNEL_SINE_TAYLOR_01, FRESNEL_SINE_TAYLOR_02,
        FRESNEL_SINE_TAYLOR_03, FRESNEL_SINE_TAYLOR_04, FRESNEL_SINE_TAYLOR_05,
        FRESNEL_SINE_TAYLOR_06, FRESNEL_SINE_TAYLOR_07, FRESNEL_SINE_TAYLOR_08,
        FRESNEL_SINE_TAYLOR_09, FRESNEL_SINE_TAYLOR_10, FRESNEL_SINE_TAYLOR_11,
        FRESNEL_SINE_TAYLOR_12, FRESNEL_SINE_TAYLOR_13, FRESNEL_SINE_TAYLOR_14,
        FRESNEL_SINE_TAYLOR_15, FRESNEL_SINE_TAYLOR_16, FRESNEL_SINE_TAYLOR_17,
        FRESNEL_SINE_TAYLOR_18, FRESNEL_SINE_TAYLOR_19, FRESNEL_SINE_TAYLOR_20,
        FRESNEL_SINE_TAYLOR_21, FRESNEL_SINE_TAYLOR_22, FRESNEL_SINE_TAYLOR_23,
        FRESNEL_SINE_TAYLOR_24, FRESNEL_SINE_TAYLOR_25, FRESNEL_SINE_TAYLOR_26
    };

    /* Variables for S(x) and powers of x, respectively. */
    double sx;
    double arg = x*x;
    double x4 = arg*arg;
    double EPS = 1.0e-8;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        int i = 0;
        double term = arg*FRESNEL_SINE_TAYLOR_COEFFICIENTS[0];
        sx = term;
        while (term > EPS){
            /* Odd terms have a negative coefficients.
               Compute two terms at a time to compare with EPS. */
            i += 1;
            arg *= x4;
            term = arg*FRESNEL_SINE_TAYLOR_COEFFICIENTS[i];
            sx += term;

            // Compute even term.
            i += 1;
            arg *= x4;
            term = arg*FRESNEL_SINE_TAYLOR_COEFFICIENTS[i];
            sx += term;
        }
        return x*sx;
    }
    else if (arg < 1.0e16) {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        cosarg *= arg;
        arg *= arg;
        sinarg *= arg;

        cosarg *= FRESNEL_SINE_ASYM_00 + arg*(
                    FRESNEL_SINE_ASYM_02 + arg*(
                        FRESNEL_SINE_ASYM_04 + FRESNEL_SINE_ASYM_06*arg
                    )
                );
        sinarg *= FRESNEL_SINE_ASYM_01 + arg*(
                    FRESNEL_SINE_ASYM_03 + arg*(
                        FRESNEL_SINE_ASYM_05 + FRESNEL_SINE_ASYM_07*arg
                    )
                );

        sx = cosarg + sinarg;
        sx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        return sx + ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    else {
        /* For large values, return the limit of S(x) as x -> +/- infinity. */
        return ((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
}

double Fresnel_Sine_Heald_Rational_EPS_Minus_Three(double x)
{
    double A, R, a, b, c, d;
    x *= SQRT_2_BY_PI;

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

    return SQRT_PI_BY_8 - R*cos(A);
}

double Fresnel_Sine_Heald_Rational_EPS_Minus_Four(double x)
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

    return sgn_x*(SQRT_PI_BY_8 - R*cos(A));
}

double Fresnel_Sine_Heald_Rational_EPS_Minus_Six(double x)
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

    return sgn_x*(SQRT_PI_BY_8 - R*cos(A));
}

double Fresnel_Sine_Heald_Rational_EPS_Minus_Eight(double x)
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

    return sgn_x*(SQRT_PI_BY_8 - R*cos(A));
}
