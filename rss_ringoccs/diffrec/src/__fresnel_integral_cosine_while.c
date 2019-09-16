/*  Various coefficients and constants defined here.                          */
#include "__math_constants.h"

/*  Various Fresnel integral functions declared here.                         */
#include "__fresnel_integrals.h"

float Fresnel_Cosine_While_to_Asymptotic_Float(float x)
{
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

double Fresnel_Cosine_While_to_Asymptotic_Double(double x)
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

long double Fresnel_Cosine_While_to_Asymptotic_Long_Long_Double(long double x)
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