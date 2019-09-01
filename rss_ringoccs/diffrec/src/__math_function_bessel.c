#include "__math_functions.h"

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