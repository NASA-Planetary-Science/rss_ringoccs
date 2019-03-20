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

#include <math.h>
#include <complex.h>

/* Define Miscellaneous Constants. */
#define SQRT_PI_BY_8 0.6266570686577501
#define SQRT_PI_BY_2 1.2533141373155001
#define SQRT_2_BY_PI 0.7978845608028654
#define PI_BY_TWO 1.5707963267948966
#define SQRT_2 1.4142135623730951

/* Define Coefficients for the Fresnel Cosine Taylor Expansion. */
#define FRESNEL_COSINE_TAYLOR_00 1.0
#define FRESNEL_COSINE_TAYLOR_01 -0.1
#define FRESNEL_COSINE_TAYLOR_02 0.004629629629629629
#define FRESNEL_COSINE_TAYLOR_03 -0.00010683760683760684
#define FRESNEL_COSINE_TAYLOR_04 1.4589169000933706e-06
#define FRESNEL_COSINE_TAYLOR_05 -1.3122532963802806e-08
#define FRESNEL_COSINE_TAYLOR_06 8.35070279514724e-11
#define FRESNEL_COSINE_TAYLOR_07 -3.9554295164585257e-13
#define FRESNEL_COSINE_TAYLOR_08 1.4483264643598138e-15
#define FRESNEL_COSINE_TAYLOR_09 -4.221407288807088e-18
#define FRESNEL_COSINE_TAYLOR_10 1.0025164934907719e-20
#define FRESNEL_COSINE_TAYLOR_11 -1.977064753877905e-23
#define FRESNEL_COSINE_TAYLOR_12 3.289260349175752e-26
#define FRESNEL_COSINE_TAYLOR_13 -4.678483515518486e-29
#define FRESNEL_COSINE_TAYLOR_14 5.754191643982172e-32
#define FRESNEL_COSINE_TAYLOR_15 -6.180307588222796e-35
#define FRESNEL_COSINE_TAYLOR_16 5.846755007468836e-38
#define FRESNEL_COSINE_TAYLOR_17 -4.908923964523423e-41
#define FRESNEL_COSINE_TAYLOR_18 3.6824935154611457e-44
#define FRESNEL_COSINE_TAYLOR_19 -2.483069097454912e-47
#define FRESNEL_COSINE_TAYLOR_20 1.513107949541217e-50
#define FRESNEL_COSINE_TAYLOR_21 -8.373419683872281e-54
#define FRESNEL_COSINE_TAYLOR_22 4.2267897541935526e-57
#define FRESNEL_COSINE_TAYLOR_23 -1.954102582324171e-60
#define FRESNEL_COSINE_TAYLOR_24 8.30461450592911e-64
#define FRESNEL_COSINE_TAYLOR_25 -3.255395462013028e-67
#define FRESNEL_COSINE_TAYLOR_26 1.1807618389115701e-70

/* Define Coefficients for the Fresnel Coine Asymptotic Expansion. */
#define FRESNEL_COSINE_ASYM_00 0.5
#define FRESNEL_COSINE_ASYM_01 -0.25
#define FRESNEL_COSINE_ASYM_02 -0.375
#define FRESNEL_COSINE_ASYM_03 0.9375
#define FRESNEL_COSINE_ASYM_04 3.281250
#define FRESNEL_COSINE_ASYM_05 -14.765625
#define FRESNEL_COSINE_ASYM_06 -81.210938
#define FRESNEL_COSINE_ASYM_07 527.87109375

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