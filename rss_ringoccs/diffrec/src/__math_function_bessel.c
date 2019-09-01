#include "__math_functions.h"

/******************************************************************************
 *  Functions:                                                                *
 *      BesselJ0_Float, BesselJ0_Double, BesselJ0_Long_Double                 *
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
 *                      _____                                                 *
 *          J_0(x)  =   \      (-1)^n x^2n /                                  *
 *                      /____             / (n)!^2 * 4^n                      *
 *                      n = 0                                                 *
 *                                                                            *
 *      For large arguments the asymptotic expansion is used. This is defined *
 *      by the following series:                                              *
 *                                                                            *
 *                      _____                                                 *
 *          J_0(x)  ~   \      cos(z) a_{2n} /    + sin(z) a_{2n+1} /         *
 *                      /____               / x^2n                 / x^{2n+1} *
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
float BesselJ0_Float(float x){

    x = fabsf(x);
    float arg = x*x;

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
    else if (arg < 1.0e32) {
        float sinarg, cosarg;

        arg = 1.0/arg;

        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;
        sinarg *= sinf(x - PI_BY_FOUR)/x;

        cosarg  = arg * BESSEL_J0_ASYM_06 + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= cosf(x - PI_BY_FOUR);

        return (cosarg + sinarg)*SQRT_2_BY_PI/sqrtf(x);

    }
    else {
        return 0.0;
    }
}

double BesselJ0_Double(double x){

    x = fabs(x);
    double arg = x*x;

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
    else if (arg < 1.0e32) {
        double sinarg, cosarg;

        arg = 1.0/arg;

        sinarg  = arg * BESSEL_J0_ASYM_07 + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;
        sinarg *= sin(x - PI_BY_FOUR)/x;

        cosarg  = arg * BESSEL_J0_ASYM_06 + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= cos(x - PI_BY_FOUR);

        return (cosarg + sinarg)*SQRT_2_BY_PI/sqrt(x);

    }
    else {
        return 0.0;
    }
}

long double BesselJ0_Long_Double(long double x){

    x = fabsl(x);
    long double arg = x*x;

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
    else if (arg < 1.0e32) {
        long double sinarg, cosarg;

        arg = 1.0/arg;

        sinarg  = arg * BESSEL_J0_ASYM_09 + BESSEL_J0_ASYM_07;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_05;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_03;
        sinarg  = arg * sinarg + BESSEL_J0_ASYM_01;
        sinarg *= sinl(x - PI_BY_FOUR)/x;

        cosarg  = arg * BESSEL_J0_ASYM_08 + BESSEL_J0_ASYM_06;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_04;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_02;
        cosarg  = arg * cosarg + BESSEL_J0_ASYM_00;
        cosarg *= cosl(x - PI_BY_FOUR);

        return (cosarg + sinarg)*SQRT_2_BY_PI/sqrtl(x);

    }
    else {
        return 0.0;
    }
}

float BesselI0_Float(float x){

    x = fabsf(x);
    float bessel_I0, arg;

    if (x < 12.0){
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
    else if (x < 87.49) {
        arg = 1.0/x;

        bessel_I0 = arg * BESSEL_I0_ASYM_04 + BESSEL_I0_ASYM_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_00;

        return bessel_I0 * expf(x)/sqrtf(TWO_PI*x);
    }
    else {
        return INFINITY;
    }
}

double BesselI0_Double(double x){

    x = fabs(x);
    double bessel_I0, arg;

    if (x < 16.0){
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
    else if (x < 709.0) {
        arg = 1.0/x;

        bessel_I0 = arg * BESSEL_I0_ASYM_06 + BESSEL_I0_ASYM_05;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_04;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_00;

        return bessel_I0 * exp(x)/sqrt(TWO_PI*x);

    }
    else {
        return INFINITY;
    }
}

long double BesselI0_Long_Double(long double x){

    x = fabsl(x);
    long double bessel_I0, arg;

    if (x < 19.0){
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
    else if (x < 11356.34) {
        arg = 1.0/x;

        bessel_I0 = arg * BESSEL_I0_ASYM_06 + BESSEL_I0_ASYM_05;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_04;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_03;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_02;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_01;
        bessel_I0 = arg * bessel_I0 + BESSEL_I0_ASYM_00;

        return bessel_I0 * expl(x)/sqrtl(TWO_PI*x);

    }
    else {
        return INFINITY;
    }
}
