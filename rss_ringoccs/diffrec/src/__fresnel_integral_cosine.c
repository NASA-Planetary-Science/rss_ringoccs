/*  Various coefficients and constants defined here.                          */
#include "__math_constants.h"

/*  Various Fresnel integral functions declared here.                         */
#include "__fresnel_integrals.h"

float Fresnel_Cosine_Taylor_to_Asymptotic_Float(float x)
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

double Fresnel_Cosine_Taylor_to_Asymptotic_Double(double x)
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

long double Fresnel_Cosine_Taylor_to_Asymptotic_Long_Double(long double x)
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
