#ifndef FRESNEL_INTEGRAL_H
#define FRESNEL_INTEGRAL_H

#include <math.h>
#include <complex.h>

/*  Fresnel Sine and Fresnel Cosine defined here.       */
#include "__fresnel_sine.h"
#include "__fresnel_cosine.h"

double complex Fresnel_Taylor_to_Asymptotic_Func(double x)
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
    double arg, arg2;
    arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        double cx, sx;
        arg2 = arg*arg;
        cx = arg2 * FRESNEL_COSINE_TAYLOR_15 + FRESNEL_COSINE_TAYLOR_14;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_13;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_12;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_11;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_10;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_09;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_08;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_07;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_06;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_05;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_04;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_03;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_02;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_01;
        cx = arg2 * cx + FRESNEL_COSINE_TAYLOR_00;

        sx = arg2 * FRESNEL_SINE_TAYLOR_15 + FRESNEL_SINE_TAYLOR_14;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_13;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_12;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_11;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_10;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_09;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_08;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_07;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_06;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_05;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_04;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_03;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_02;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_01;
        sx = arg2 * sx + FRESNEL_SINE_TAYLOR_00;

        return x*(cx + arg*sx*_Complex_I);
    }
    else if (arg < 1.0e16) {
        double complex sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        arg2 = arg*arg;
        cosarg *= arg;
        sinarg *= arg;

        sinarg *= (
            FRESNEL_COSINE_ASYM_00 + arg2*(
                FRESNEL_COSINE_ASYM_02 + arg2*(
                    FRESNEL_COSINE_ASYM_04 + FRESNEL_COSINE_ASYM_06*arg2
                )
            )
        ) + _Complex_I * arg *(
            FRESNEL_SINE_ASYM_01 + arg2*(
                FRESNEL_SINE_ASYM_03 + arg2*(
                    FRESNEL_SINE_ASYM_05 + FRESNEL_SINE_ASYM_07*arg2
                )
            )
        );

        cosarg *= arg * (
            FRESNEL_COSINE_ASYM_01 + arg2*(
                FRESNEL_COSINE_ASYM_03 + arg2*(
                    FRESNEL_COSINE_ASYM_05 + FRESNEL_COSINE_ASYM_07*arg2
                )
            )
        ) + _Complex_I *(
            FRESNEL_SINE_ASYM_00 + arg2*(
                FRESNEL_SINE_ASYM_02 + arg2*(
                    FRESNEL_SINE_ASYM_04 + FRESNEL_SINE_ASYM_06*arg2
                )
            )
        );

        return x*(cosarg + sinarg) + 
                (1.0+_Complex_I)*((x > 0) - (x < 0))*SQRT_PI_BY_8;
    }
    else {
        /* For large values, return the limit of S(x) as x -> +/- infinity. */
        return ((x > 0) - (x < 0))*SQRT_PI_BY_8*(1.0+1.0*_Complex_I);
    }
}

double complex Fresnel_While_to_Asymptotic_Func(double x)
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

double complex Fresnel_Heald_Rational_EPS_Minus_Three(double x)
{
    double A, R, a, b, c, d, sgn_x, sx, cx;
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

    cx = sgn_x*(SQRT_PI_BY_8 - R*sin(A));
    sx = sgn_x*(SQRT_PI_BY_8 - R*cos(A));

    return cx + _Complex_I*sx;
}

double complex Fresnel_Heald_Rational_EPS_Minus_Four(double x)
{
    double A, R, a, b, c, d, sgn_x, sx, cx;
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

    cx = sgn_x*(SQRT_PI_BY_8 - R*sin(A));
    sx = sgn_x*(SQRT_PI_BY_8 - R*cos(A));

    return cx + _Complex_I*sx;
}

double complex Fresnel_Heald_Rational_EPS_Minus_Six(double x)
{
    double A, R, a, b, c, d, sgn_x, cx, sx;
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

    cx = sgn_x*(SQRT_PI_BY_8 - R*sin(A));
    sx = sgn_x*(SQRT_PI_BY_8 - R*cos(A));

    return cx + _Complex_I*sx;
}

double complex Fresnel_Heald_Rational_EPS_Minus_Eight(double x)
{
    double A, R, a, b, c, d, sgn_x, cx, sx;
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

    cx = sgn_x*(SQRT_PI_BY_8 - R*sin(A));
    sx = sgn_x*(SQRT_PI_BY_8 - R*cos(A));

    return cx + _Complex_I*sx;
}

double complex Square_Well_Diffraction_Solution(double x, double a,
                                                double b, double F)
{
    double arg1 = SQRT_PI_BY_2*(a-x)/F;
    double arg2 = SQRT_PI_BY_2*(b-x)/F;
    double complex result =   Fresnel_Heald_Rational_EPS_Minus_Eight(arg2)
                            - Fresnel_Heald_Rational_EPS_Minus_Eight(arg1);
    result *= SQRT_2_BY_PI;

    return 1.0 - (0.5 - 0.5*_Complex_I)*result;
}

double complex Inverted_Square_Well_Diffraction_Solution(double x, double a,
                                                         double b, double F)
{
    double arg1 = SQRT_PI_BY_2*(a-x)/F;
    double arg2 = SQRT_PI_BY_2*(b-x)/F;
    double complex result =   Fresnel_Heald_Rational_EPS_Minus_Eight(arg2)
                            - Fresnel_Heald_Rational_EPS_Minus_Eight(arg1);
    result *= SQRT_2_BY_PI;
    result *= SQRT_2_BY_PI;

    return (0.5 - 0.5*_Complex_I)*result;
}

#endif
