#include <math.h>

/*  Various coefficients and constants defined here.                          */
#include "__math_constants.h"

/*  Various Fresnel integral functions declared here.                         */
#include "__fresnel_integrals.h"

/*----------------------Single Precision Functions----------------------------*/

float Fresnel_Cosine_Taylor_to_Asymptotic_Float(float x){

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

float Fresnel_Cosine_While_to_Asymptotic_Float(float x){
    float FRESNEL_COSINE_TAYLOR_COEFFICIENTS[30] = {
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
        FRESNEL_COSINE_TAYLOR_26, FRESNEL_COSINE_TAYLOR_27,
        FRESNEL_COSINE_TAYLOR_28, FRESNEL_COSINE_TAYLOR_29
    };

    /* Variables for S(x) and powers of x, respectively. */
    float cx;
    float arg = x*x;
    float x4 = arg*arg;
    float EPS = 1.0e-8;

    if (arg < 9.0){
        int i = 0;
        float term = arg*FRESNEL_COSINE_TAYLOR_COEFFICIENTS[0];
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
        float sinarg, cosarg;
        cosarg = cosf(arg);
        sinarg = sinf(arg);
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

float Fresnel_Cosine_Heald_Rational_EPS_Minus_Three_Float(float x){
    float A, R, a, b, c, d, sgn_x;
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinf(A));
}

float Fresnel_Cosine_Heald_Rational_EPS_Minus_Four_Float(float x){
    float A, R, a, b, c, d, sgn_x;
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinf(A));
}

float Fresnel_Cosine_Heald_Rational_EPS_Minus_Six_Float(float x){
    float A, R, a, b, c, d, sgn_x;
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinf(A));
}

float Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Float(float x){
    float A, R, a, b, c, d, sgn_x;
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinf(A));
}

/*----------------------Double Precision Functions----------------------------*/

double Fresnel_Cosine_Taylor_to_Asymptotic_Double(double x){

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

double Fresnel_Cosine_While_to_Asymptotic_Double(double x){
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

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Three_Double(double x){
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

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Four_Double(double x){
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

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Six_Double(double x){
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

double Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Double(double x){
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

/*--------------------Long Double Precision Functions-------------------------*/

long double Fresnel_Cosine_Taylor_to_Asymptotic_Long_Double(long double x){

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

long double Fresnel_Cosine_While_to_Asymptotic_Long_Long_Double(long double x){
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
    long double cx;
    long double arg = x*x;
    long double x4 = arg*arg;
    long double EPS = 1.0e-8;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        int i = 0;
        long double term = arg*FRESNEL_COSINE_TAYLOR_COEFFICIENTS[0];
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
        long double sinarg, cosarg;
        cosarg = cosl(arg);
        sinarg = sinl(arg);
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

long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Three_Long_Double(
    long double x
){
    long double A, R, a, b, c, d, sgn_x;
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinl(A));
}

long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Four_Long_Double(
    long double x
){
    long double A, R, a, b, c, d, sgn_x;
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinl(A));
}

long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Six_Long_Double(
    long double x
){
    long double A, R, a, b, c, d, sgn_x;
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinl(A));
}

long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Long_Double(
    long double x
){
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

    return sgn_x*(SQRT_PI_BY_8 - R*sinl(A));
}