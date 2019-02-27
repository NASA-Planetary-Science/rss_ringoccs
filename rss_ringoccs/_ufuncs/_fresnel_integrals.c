/*******************************************************************************
 *                             Fresnel Integrals                               *
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
 *                                                                             *
 *           x                                                                 *
 *           -                                                                 *
 *          | |                                                                *
 * C(x) =   |   cos(t^2) dt                                                    *
 *        | |                                                                  *
 *         -                                                                   *
 *         0                                                                   *
 *******************************************************************************
 * It is very common for a pi/2 to be placed inside the sine and cosine terms, *
 * and thus in translating one would need to scale x by sqrt(2/pi) and scale   *
 * the results by sqrt(pi/2). Several of the algorithms implemented compute    *
 * the latter definitions. We have taken this into account in our function and *
 * return the values corresponding to the equations above.                     *
 *******************************************************************************
 *                              DEFINED FUNCTIONS                              *
 *******************************************************************************
 *  Fresnel_Sine_Taylor_to_Asymptotic:                                         *
 *      Description.                                                           *
 *  Fresnel_Sine_While_to_Asymptocic                                           *
 *      Description.                                                           *
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Three                                *
 *      Description.                                                           *
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Four                                 *
 *      Description.                                                           *
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Six                                  *
 *      Description.                                                           *
 *  Fresnel_Sine_Heald_Rational_EPS_Minus_Eight                                *
 *      Description.                                                           *
 *******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                                *
 *  Date:       Febuary 26, 2019                                               *
 ******************************************************************************/

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <math.h>
#include <complex.h>
#include "../../include/Python.h"
#include "../../include/ndarraytypes.h"
#include "../../include/ufuncobject.h"

/* Define Miscellaneous Constants. */
#define SQRT_PI_BY_8 0.6266570686577501
#define SQRT_PI_BY_2 1.2533141373155001
#define SQRT_2_BY_PI 0.7978845608028654
#define PI_BY_TWO 1.5707963267948966

/*
 *  DEFINE COEFFICIENTS FOR FRESNEL SINE INTEGRAL.
 */

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

/* Define Coefficients for the Fresnel Sine Asymptotic Expansion. */
#define FRESNEL_SINE_ASYM_00 -0.5
#define FRESNEL_SINE_ASYM_01 -0.25
#define FRESNEL_SINE_ASYM_02 0.375
#define FRESNEL_SINE_ASYM_03 0.9375

/*
 *  DEFINE COEFFICIENTS FOR FRESNEL COSINE INTEGRAL.
 */

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
#define FRESNEL_COSINE_TAYLOR_10 1.0025164934907719e-20f
#define FRESNEL_COSINE_TAYLOR_11 -1.977064753877905e-23
#define FRESNEL_COSINE_TAYLOR_12 3.289260349175752e-26
#define FRESNEL_COSINE_TAYLOR_13 -4.6784835155184856e-29

/* Define Coefficients for the Fresnel Coine Asymptotic Expansion. */
#define FRESNEL_COSINE_ASYM_00 0.5
#define FRESNEL_COSINE_ASYM_01 -0.25
#define FRESNEL_COSINE_ASYM_02 -0.375
#define FRESNEL_COSINE_ASYM_03 0.9375

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
#define FRESNEL_HEALD_RATIONAL_EPS_8_A09 1.0000000
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

static PyMethodDef _fresnel_integrals_methods[] = {{NULL, NULL, 0, NULL}};
/*-----------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.  *
 * The are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.               *
 *----------------------------------------------------------------------------*/
double Fresnel_Sine_Taylor_to_Asymptotic_Func(double x)
{
    /* Variables for S(x) and powers of x, respectively. */
    double sx;
    double arg = x*x;

    /* For small x use the Taylor expansion to compute C(x). For larger x,  *
     * use the asymptotic expansion. For values near 3.076, accuracy of 5   *
     * decimals is guaranteed. Higher precicion outside this region.        */
    if (arg < 9.0){
        double arg_sq = arg*arg;
        if (arg < 1.0){
            sx = arg_sq * FRESNEL_SINE_TAYLOR_03 + FRESNEL_SINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_00;
            sx *= arg;
            return sx*x;
        }
        else if (arg < 4.0){
            sx = arg_sq * FRESNEL_SINE_TAYLOR_07 + FRESNEL_SINE_TAYLOR_06;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_05;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_04;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_03;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_00;
            sx *= arg;
            return sx*x;
        }
        else{
            sx = arg_sq * FRESNEL_SINE_TAYLOR_13 + FRESNEL_SINE_TAYLOR_12;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_11;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_10;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_09;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_08;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_07;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_06;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_05;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_04;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_03;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_02;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_01;
            sx = arg_sq * sx + FRESNEL_SINE_TAYLOR_00;
            sx *= arg;
            return sx*x;
        }
    }
    else {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        cosarg *= arg;
        arg *= arg;
        sinarg *= arg;
        cosarg *= FRESNEL_SINE_ASYM_02*arg + FRESNEL_SINE_ASYM_00;
        sinarg *= FRESNEL_SINE_ASYM_03*arg + FRESNEL_SINE_ASYM_01;

        sx = cosarg + sinarg;
        sx *= x;
        if (x > 0){
            return sx+SQRT_PI_BY_8;
        }
        else {
            return sx-SQRT_PI_BY_8;
        }
    }
}

double Fresnel_Sine_While_to_Asymptotic_Func(double x)
{
    double FRESNEL_SINE_TAYLOR_COEFFICIENTS[14] = {
        FRESNEL_SINE_TAYLOR_00, FRESNEL_SINE_TAYLOR_01, FRESNEL_SINE_TAYLOR_02,
        FRESNEL_SINE_TAYLOR_03, FRESNEL_SINE_TAYLOR_04, FRESNEL_SINE_TAYLOR_05,
        FRESNEL_SINE_TAYLOR_06, FRESNEL_SINE_TAYLOR_07, FRESNEL_SINE_TAYLOR_08,
        FRESNEL_SINE_TAYLOR_09, FRESNEL_SINE_TAYLOR_10, FRESNEL_SINE_TAYLOR_11,
        FRESNEL_SINE_TAYLOR_12, FRESNEL_SINE_TAYLOR_13
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
    else {
        double sinarg, cosarg;
        cosarg = cos(arg);
        sinarg = sin(arg);
        arg = 1.0/arg;
        cosarg *= arg;
        arg *= arg;
        sinarg *= arg;
        cosarg *= FRESNEL_SINE_ASYM_02*arg + FRESNEL_SINE_ASYM_00;
        sinarg *= FRESNEL_SINE_ASYM_03*arg + FRESNEL_SINE_ASYM_01;

        sx = cosarg + sinarg;
        sx *= x;
        if (x > 0){
            return sx+SQRT_PI_BY_8;
        }
        else {
            return sx-SQRT_PI_BY_8;
        }
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
    double A, R, a, b, c, d;
    x *= SQRT_2_BY_PI;

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

    return SQRT_PI_BY_8 - R*cos(A);
}

double Fresnel_Sine_Heald_Rational_EPS_Minus_Six(double x)
{
    double A, R, a, b, c, d;
    x *= SQRT_2_BY_PI;

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

    return SQRT_PI_BY_8 - R*cos(A);
}

/*---------------------------DEFINE PYTHON FUNCTIONS--------------------------*
 * This contains the Numpy-C and Python-C API parts that allow for the above  *
 * functions to be called in Python. Numpy arrays, as well as floating point  *
 * and integer valued arguments may then be passed into these functions for   *
 * improvement in performance, as opposed to the routines written purely in   *
 * Python. Successful compiling requires the Numpy and Python header files.   *
 *----------------------------------------------------------------------------*/
static void double_fresnelsin_taylor_to_asymp(char **args, npy_intp *dimensions,
                                              npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in1 = args[0];
    char *out1 = args[1];
    npy_intp in1_step = steps[0];
    npy_intp out1_step = steps[1];

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        *((double *)out1) = Fresnel_Sine_Taylor_to_Asymptotic_Func(
            *(double *)in1
        );
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

static void double_fresnelsin_while_to_asymp(char **args, npy_intp *dimensions,
                                             npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in1 = args[0];
    char *out1 = args[1];
    npy_intp in1_step = steps[0];
    npy_intp out1_step = steps[1];

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        *((double *)out1) = Fresnel_Sine_While_to_Asymptotic_Func(
            *(double *)in1
        );
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

static void double_fresnel_heald_eps_3(char **args, npy_intp *dimensions,
                                       npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in1 = args[0];
    char *out1 = args[1];
    npy_intp in1_step = steps[0];
    npy_intp out1_step = steps[1];

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        *((double *)out1) = Fresnel_Sine_Heald_Rational_EPS_Minus_Three(
            *(double *)in1
        );
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

static void double_fresnel_heald_eps_4(char **args, npy_intp *dimensions,
                                       npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in1 = args[0];
    char *out1 = args[1];
    npy_intp in1_step = steps[0];
    npy_intp out1_step = steps[1];

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        *((double *)out1) = Fresnel_Sine_Heald_Rational_EPS_Minus_Four(
            *(double *)in1
        );
        /*END main ufunc computation*/

        in1 += in1_step;
        out1 += out1_step;
    }
}

/* Define pointers to the C functions. */
PyUFuncGenericFunction fresnel_sin_1[1] = {&double_fresnelsin_taylor_to_asymp};
PyUFuncGenericFunction fresnel_sin_2[1] = {&double_fresnelsin_while_to_asymp};
PyUFuncGenericFunction fresnel_sin_3[1] = {&double_fresnel_heald_eps_3};
PyUFuncGenericFunction fresnel_sin_4[1] = {&double_fresnel_heald_eps_4};

/* Input and return types for double input and out. */
static char double_double_types[2] = {NPY_DOUBLE, NPY_DOUBLE};
static void *PyuFunc_data[1] = {NULL};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_fresnel_integrals",
    NULL,
    -1,
    _fresnel_integrals_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__fresnel_integrals(void)
{
    PyObject *fresnel_sin_taylor_to_asymptotic;
    PyObject *fresnel_sin_while_to_asymptotic;
    PyObject *fresnel_sin_heald_eps_three;
    PyObject *fresnel_sin_heald_eps_four;
    PyObject *fresnel_sin_5;

    PyObject *m, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }

    import_array();
    import_umath();

    fresnel_sin_taylor_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_sin_1, PyuFunc_data, double_double_types,
        1, 1, 1, PyUFunc_None,
        "fresnel_sin_taylor_to_asymptotic",
        "fresnel_sin_1_docstring", 0
    );

    fresnel_sin_while_to_asymptotic = PyUFunc_FromFuncAndData(
        fresnel_sin_2, PyuFunc_data, double_double_types,
        1, 1, 1, PyUFunc_None, 
        "fresnel_sin_while_to_asymptotic",
        "fresnel_sin_while_to_asymptotic_docstring", 0
    );

    fresnel_sin_heald_eps_three = PyUFunc_FromFuncAndData(
        fresnel_sin_3, PyuFunc_data, double_double_types,
        1, 1, 1, PyUFunc_None, 
        "fresnel_sin_heald_eps_three",
        "fresnel_sin_heald_eps_three_docstring", 0
    );
    fresnel_sin_heald_eps_four = PyUFunc_FromFuncAndData(
        fresnel_sin_4, PyuFunc_data, double_double_types,
        1, 1, 1, PyUFunc_None, 
        "fresnel_sin_heald_eps_three",
        "fresnel_sin_heald_eps_three_docstring", 0
    );

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin_taylor_to_asymptotic",
                         fresnel_sin_taylor_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_sin_while_to_asymptotic",
                         fresnel_sin_while_to_asymptotic);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_three",
                         fresnel_sin_heald_eps_three);
    PyDict_SetItemString(d, "fresnel_sin_heald_eps_four",
                         fresnel_sin_heald_eps_four);
    Py_DECREF(fresnel_sin_taylor_to_asymptotic);
    Py_DECREF(fresnel_sin_while_to_asymptotic);
    Py_DECREF(fresnel_sin_heald_eps_three);
    Py_DECREF(fresnel_sin_heald_eps_four);
    return m;
}
#else
PyMODINIT_FUNC init__funcs(void)
{
    PyObject *fresnel_sin_1;
    PyObject *fresnel_sin_2;
    PyObject *fresnel_sin_3;
    PyObject *fresnel_sin_4;
    PyObject *fresnel_sin_5;

    PyObject *m, *d;

    m = Py_InitModule("__funcs", _fresnel_integrals_methods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    fresnel_sin_1 = PyUFunc_FromFuncAndData(fresnel_sin_funcs_1, PyuFunc_data,
                                            double_double_types, 1, 1, 1,
                                            PyUFunc_None, "fresnel_sin_1",
                                            "fresnel_sin_1_docstring", 0);

    fresnel_sin_2 = PyUFunc_FromFuncAndData(fresnel_sin_funcs_2, PyuFunc_data,
                                            double_double_types, 1, 1, 1,
                                            PyUFunc_None, "fresnel_sin_2",
                                            "fresnel_sin_2_docstring", 0);

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "fresnel_sin_1", fresnel_sin_1);
    PyDict_SetItemString(d, "fresnel_sin_2", fresnel_sin_2);
    Py_DECREF(fresnel_sin_1);
    Py_DECREF(fresnel_sin_2);
}
#endif