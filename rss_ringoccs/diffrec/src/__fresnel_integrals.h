/*******************************************************************************
 *                           Fresnel Integrals                                 *
 *******************************************************************************
 *  This C program contains several algorithms for the computation of the      *
 *  Fresnel Cosine integral. This is meant to test the accuracy and efficiency *
 *  of the various algorithms. This file is not included in the setup.py file, *
 *  and thus will not be compiled. The fresnel sine and cosine functions that  *
 *  are used in the _special_functions.c file are the functions that had the   *
 *  best combination of speed and accuracy. This file is kept so that users    *
 *  may experiment with the various known algorithms.                          *
 *******************************************************************************
 *  We define the Fresnel Cosine Integrals as follows:                         *
 *                 x                                                           *
 *                 -                                                           *
 *                | |                                                          *
 *       C(x) =   |   cos(t^2) dt                                              *
 *              | |                                                            *
 *               -                                                             *
 *               0                                                             *
 *******************************************************************************
 *  It is very common for a pi/2 to be placed inside the cosine term,          *
 *  and thus in translating one would need to scale x by sqrt(2/pi) and scale  *
 *  the results by sqrt(pi/2). Several of the algorithms implemented compute   *
 *  the latter definition. We have taken this into account in our functions    *
 *  and return the values corresponding to the equations above.                *
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
 *  Date:       August 30, 2019                                                *
 ******************************************************************************/
#ifndef RSS_RINGOCCS_FRESNEL_INTEGRALS_H
#define RSS_RINGOCCS_FRESNEL_INTEGRALS_H

#include <math.h>
#include <complex.h>

/*  Various coefficients and constants defined here.                          */
#include "__math_constants.h"

/******************************************************************************
 *----------------------------Fresnel Cosine----------------------------------*
 ******************************************************************************/

/*----------------------Single Precision Functions----------------------------*/
extern float Fresnel_Sine_Taylor_to_Asymptotic_Float(float x);

extern float Fresnel_Sine_While_to_Asymptotic_Float(float x);

extern float Fresnel_Sine_Heald_Rational_EPS_Minus_Three_Float(float x);

extern float Fresnel_Sine_Heald_Rational_EPS_Minus_Four_Float(float x);

extern float Fresnel_Sine_Heald_Rational_EPS_Minus_Six_Float(float x);

extern float Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Float(float x);

/*----------------------Double Precision Functions----------------------------*/
extern double Fresnel_Sine_Taylor_to_Asymptotic_Double(double x);

extern double Fresnel_Sine_While_to_Asymptotic_Double(double x);

extern double Fresnel_Sine_Heald_Rational_EPS_Minus_Three_Double(double x);

extern double Fresnel_Sine_Heald_Rational_EPS_Minus_Four_Double(double x);

extern double Fresnel_Sine_Heald_Rational_EPS_Minus_Six_Double(double x);

extern double Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Double(double x);

/*--------------------Long Double Precision Functions-------------------------*/
extern long double Fresnel_Sine_Taylor_to_Asymptotic_Long_Double(
    long double x
);

extern long double Fresnel_Sine_While_to_Asymptotic_Long_Long_Double(
    long double x
);

extern long double Fresnel_Sine_Heald_Rational_EPS_Minus_Three_Long_Double(
    long double x
);

extern long double Fresnel_Sine_Heald_Rational_EPS_Minus_Four_Long_Double(
    long double x
);

extern long double Fresnel_Sine_Heald_Rational_EPS_Minus_Six_Long_Double(
    long double x
);

extern long double Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Long_Double(
    long double x
);

/******************************************************************************
 *-----------------------------Fresnel Sine-----------------------------------*
 ******************************************************************************/

/*----------------------Single Precision Functions----------------------------*/
extern float Fresnel_Cosine_Taylor_to_Asymptotic_Float(float x);

extern float Fresnel_Cosine_While_to_Asymptotic_Float(float x);

extern float Fresnel_Cosine_Heald_Rational_EPS_Minus_Three_Float(float x);

extern float Fresnel_Cosine_Heald_Rational_EPS_Minus_Four_Float(float x);

extern float Fresnel_Cosine_Heald_Rational_EPS_Minus_Six_Float(float x);

extern float Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Float(float x);

/*----------------------Double Precision Functions----------------------------*/
extern double Fresnel_Cosine_Taylor_to_Asymptotic_Double(double x);

extern double Fresnel_Cosine_While_to_Asymptotic_Double(double x);

extern double Fresnel_Cosine_Heald_Rational_EPS_Minus_Three_Double(double x);

extern double Fresnel_Cosine_Heald_Rational_EPS_Minus_Four_Double(double x);

extern double Fresnel_Cosine_Heald_Rational_EPS_Minus_Six_Double(double x);

extern double Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Double(double x);

/*--------------------Long Double Precision Functions-------------------------*/
extern long double Fresnel_Cosine_Taylor_to_Asymptotic_Long_Double(
    long double x
);

extern long double Fresnel_Cosine_While_to_Asymptotic_Long_Long_Double(
    long double x
);

extern long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Three_Long_Double(
    long double x
);

extern long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Four_Long_Double(
    long double x
);

extern long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Six_Long_Double(
    long double x
);

extern long double Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Long_Double(
    long double x
);

/******************************************************************************
 *------------------------Complex Fresnel Integral----------------------------*
 ******************************************************************************/

extern complex double Fresnel_Taylor_to_Asymptotic_Double(double x);

extern complex double Fresnel_While_to_Asymptotic_Func(double x);

extern complex double Fresnel_Heald_Rational_EPS_Minus_Three_Func(double x);

extern complex double Fresnel_Heald_Rational_EPS_Minus_Four_Func(double x);

extern complex double Fresnel_Heald_Rational_EPS_Minus_Six_Func(double x);

extern complex double Fresnel_Heald_Rational_EPS_Minus_Eight_Func(double x);

#endif
