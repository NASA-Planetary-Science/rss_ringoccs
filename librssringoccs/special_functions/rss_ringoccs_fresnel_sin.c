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
 *                  ___                                                      *
 *                  \      (-1)^n x^(4n+3)/                                    *
 *        S(x)=     /__                /  (4n+3)(2n+1)!                      *
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
 *                             ___                                           *
 *                             \                                               *
 *        S(x) = sqrt(pi/i) -  /__  (-1)^n (a_n(x)sin(x^2) + b_n(x)cos(x^2)) *
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

/*  The C Standard Library header for math functions and more found here.     */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Complex variables and functions defined here.                             */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/* Define Coefficients for the Fresnel Sine Taylor Expansion.                 */
#define FRESNEL_SINE_TAYLOR_00  0.33333333333333333333333333333
#define FRESNEL_SINE_TAYLOR_01 -2.38095238095238095238095238095e-2
#define FRESNEL_SINE_TAYLOR_02  7.57575757575757575757575757576e-4
#define FRESNEL_SINE_TAYLOR_03 -1.32275132275132275132275132275e-5
#define FRESNEL_SINE_TAYLOR_04  1.45038522231504687645038522232e-7
#define FRESNEL_SINE_TAYLOR_05 -1.08922210371485733804574384285e-9
#define FRESNEL_SINE_TAYLOR_06  5.94779401363763503681199154450e-12
#define FRESNEL_SINE_TAYLOR_07 -2.46682701026445692771004257606e-14
#define FRESNEL_SINE_TAYLOR_08  8.03273501241577360913984452289e-17
#define FRESNEL_SINE_TAYLOR_09 -2.10785519144213582486050800945e-19
#define FRESNEL_SINE_TAYLOR_10  4.55184675892820028624362194733e-22
#define FRESNEL_SINE_TAYLOR_11 -8.23014929921422135684449347133e-25
#define FRESNEL_SINE_TAYLOR_12  1.26410789889891635219506925867e-27
#define FRESNEL_SINE_TAYLOR_13 -1.66976179341737202698649397027e-30
#define FRESNEL_SINE_TAYLOR_14  1.91694286210978253077267196219e-33
#define FRESNEL_SINE_TAYLOR_15 -1.93035720881510785655551537411e-36
#define FRESNEL_SINE_TAYLOR_16  1.71885606280178362396819126766e-39
#define FRESNEL_SINE_TAYLOR_17 -1.36304126177913957635067836351e-42
#define FRESNEL_SINE_TAYLOR_18  9.68728023887076175384366004096e-46
#define FRESNEL_SINE_TAYLOR_19 -6.20565791963739670594197460729e-49
#define FRESNEL_SINE_TAYLOR_20  3.60157930981012591661339989697e-52
#define FRESNEL_SINE_TAYLOR_21 -1.90254122728987952723942026864e-55
#define FRESNEL_SINE_TAYLOR_22  9.18642950239868569596123672835e-59
#define FRESNEL_SINE_TAYLOR_23 -4.07013527785325672297810283986e-62
#define FRESNEL_SINE_TAYLOR_24  1.66058051345108993284425792700e-65
#define FRESNEL_SINE_TAYLOR_25 -6.25918411694871134024677459636e-69
#define FRESNEL_SINE_TAYLOR_26  2.18621042295388572102809768805e-72
#define FRESNEL_SINE_TAYLOR_27 -7.09571739181805357327043566663e-76
#define FRESNEL_SINE_TAYLOR_28  2.14564844309633852738645079818e-79
#define FRESNEL_SINE_TAYLOR_29 -6.05939744697137480782877578571e-83
#define FRESNEL_SINE_TAYLOR_30  1.60173329821314496897157652161e-86

/* Define Coefficients for the Fresnel Sine Asymptotic Expansion.             */
#define FRESNEL_SINE_ASYM_00    -0.50
#define FRESNEL_SINE_ASYM_01    -0.250
#define FRESNEL_SINE_ASYM_02     0.3750
#define FRESNEL_SINE_ASYM_03     0.93750
#define FRESNEL_SINE_ASYM_04    -3.281250
#define FRESNEL_SINE_ASYM_05    -14.7656250
#define FRESNEL_SINE_ASYM_06     81.21093750
#define FRESNEL_SINE_ASYM_07     527.871093750
#define FRESNEL_SINE_ASYM_08    -3959.0332031250
#define FRESNEL_SINE_ASYM_09    -33651.78222656250
#define FRESNEL_SINE_ASYM_10     319691.931152343750
#define FRESNEL_SINE_ASYM_11     3356765.2770996093750
#define FRESNEL_SINE_ASYM_12    -38602800.68664550781250
#define FRESNEL_SINE_ASYM_13    -482535008.583068847656250
#define FRESNEL_SINE_ASYM_14     6514222615.8714294433593750

float rssringoccs_Float_Fresnel_Sin(float x)
{
    /* Variables for S(x) and powers of x, respectively. */
    float sx;
    float arg = x*x;
    float sinarg, cosarg, cos_x_squared, sin_x_squared;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 9.0)
    {
        x *= arg;
        arg *= arg;
        sx = arg * FRESNEL_SINE_TAYLOR_16 + FRESNEL_SINE_TAYLOR_15;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_14;
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
        sx = sx*x;
    }
    else if (arg < 1.0e16)
    {
        cos_x_squared = rssringoccs_Float_Cos(arg);
        sin_x_squared = rssringoccs_Float_Sin(arg);

        arg = 1.0/arg;
        cos_x_squared *= arg;
        arg *= arg;
        sin_x_squared *= arg;

        cosarg  = arg * FRESNEL_SINE_ASYM_08 + FRESNEL_SINE_ASYM_06;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_04;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_02;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_00;
        cosarg *= cos_x_squared;

        sinarg  = arg * FRESNEL_SINE_ASYM_09 + FRESNEL_SINE_ASYM_07;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_05;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_03;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_01;
        sinarg *= sin_x_squared;

        sx = cosarg + sinarg;
        sx *= x;

        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        sx = sx + ((x > 0) - (x < 0))*rssringoccs_Sqrt_Pi_By_Eight;
    }

    /*  For large values, return the limit of S(x) as x -> +/- infinity.      */
    else
        sx = ((x > 0) - (x < 0))*rssringoccs_Sqrt_Pi_By_Eight;

    return sx;
}

double rssringoccs_Double_Fresnel_Sin(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double sx;
    double sinarg, cosarg, cos_x_squared, sin_x_squared;
    double arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 11.68){
        x *= arg;
        arg *= arg;
        sx = arg * FRESNEL_SINE_TAYLOR_22 + FRESNEL_SINE_TAYLOR_21;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_20;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_19;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_18;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_17;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_16;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_15;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_14;
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
        sx = sx*x;
    }
    else if (arg < 1.0e16)
    {
        cos_x_squared = rssringoccs_Double_Cos(arg);
        sin_x_squared = rssringoccs_Double_Sin(arg);

        arg = 1.0/arg;
        cos_x_squared *= arg;
        arg *= arg;
        sin_x_squared *= arg;

        cosarg  = arg * FRESNEL_SINE_ASYM_08 + FRESNEL_SINE_ASYM_06;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_04;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_02;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_00;
        cosarg *= cos_x_squared;

        sinarg  = arg * FRESNEL_SINE_ASYM_09 + FRESNEL_SINE_ASYM_07;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_05;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_03;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_01;
        sinarg *= sin_x_squared;

        sx = cosarg + sinarg;
        sx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        sx = sx + ((x > 0) - (x < 0))*rssringoccs_Sqrt_Pi_By_Eight;
    }
    /*  For large values, return the limit of S(x) as x -> +/- infinity.      */
    else
        sx = ((x > 0) - (x < 0))*rssringoccs_Sqrt_Pi_By_Eight;

    return sx;
}

long double rssringoccs_LDouble_Fresnel_Sin(long double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    long double sx;
    long double arg = x*x;
    long double sinarg, cosarg, cos_x_squared, sin_x_squared;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region. When   *
     * |x| > 1.e8, S(x) returns +/- sqrt(pi/8) to 8 decimals.               */
    if (arg < 11.68)
    {
        x *= arg;
        arg *= arg;
        sx = arg * FRESNEL_SINE_TAYLOR_24 + FRESNEL_SINE_TAYLOR_23;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_22;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_21;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_20;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_19;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_18;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_17;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_16;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_15;
        sx = arg * sx + FRESNEL_SINE_TAYLOR_14;
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
        sx = sx*x;
    }
    else if (arg < 1.0e16)
    {
        cos_x_squared = rssringoccs_LDouble_Cos(arg);
        sin_x_squared = rssringoccs_LDouble_Sin(arg);

        arg = 1.0/arg;
        cos_x_squared *= arg;
        arg *= arg;
        sin_x_squared *= arg;

        cosarg  = arg * FRESNEL_SINE_ASYM_08 + FRESNEL_SINE_ASYM_06;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_04;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_02;
        cosarg  = arg * cosarg + FRESNEL_SINE_ASYM_00;
        cosarg *= cos_x_squared;

        sinarg  = arg * FRESNEL_SINE_ASYM_09 + FRESNEL_SINE_ASYM_07;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_05;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_03;
        sinarg  = arg * sinarg + FRESNEL_SINE_ASYM_01;
        sinarg *= sin_x_squared;

        sx = cosarg + sinarg;
        sx *= x;
        /*  (x > 0) - (x < 0) is a quick way to return sign(x) and avoids an  *
         *  expensive if-then statement. Output for the asymptotic expansion  *
         *  is f(|x|) + sign(x) * sqrt(pi/8). Error goes like 1/x^15.         */
        sx = sx + ((x > 0) - (x < 0))*rssringoccs_Sqrt_Pi_By_Eight;
    }

    /* For large values, return the limit of S(x) as x -> +/- infinity.       */
    else
        sx = ((x > 0) - (x < 0))*rssringoccs_Sqrt_Pi_By_Eight;

    return sx;
}
