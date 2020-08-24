/******************************************************************************
 *                   Double Slit Fraunhofer Diffraction                       *
 ******************************************************************************
 *  These functions are used for modeling diffraction patterns from a double  * 
 *  slit when the Fraunhofer approximation is justified.                      *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  Double_Slit_Fraunhofer_Diffraction                                        *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from two slits.                       *
 *  Arguments:                                                                *
 *      x (float, double, or long double):                                    *
 *          A real number, the location on the x-axis of the screen.          *
 *      z (float, double, or long double):                                    *
 *          A real number, the location on the x-axis of the screen.          *
 *  Output:                                                                   *
 *      bessel_J0:                                                            *
 *          The Bessel function J_0(x).                                       *
 *  Method:                                                                   *
 *      For small values, the Taylor expansion is used. The J_0(x) function   *
 *      can be defined by the following series:                               *
 *                                                                            *
 *                      ___                                                 *
 *                      \      (-1)^n x^2n /                                  *
 *         J_0(x)  =    /__             / (n)!^2 * 4^n                      *
 *                      n = 0                                                 *
 *                                                                            *
 *      For large arguments the asymptotic expansion is used. This is defined *
 *      by the following series:                                              *
 *                                                                            *
 *                      ___                                                 *
 *                      \      cos(z) a_{2n} /    + sin(z) a_{2n+1} /         *
 *          J_0(x)  ~   /__               / x^2n                 / x^{2n+1} *
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
 ******************************************************************************/

#include "special_functions.h"
#include "_math_constants.h"
#include <math.h>

float Double_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a,
                                               float d, float lambda)
{
    float sin_theta = x/sqrtf(x*x + z*z);
    a = a/lambda;
    d = d/lambda;

    float var_1  = Sinc_Float(a*sin_theta);
    var_1       *= var_1;

    float var_2  = sinf(TWO_PI*d*sin_theta);
    var_2       *= var_2;

    float var_3  = 2.0*sinf(ONE_PI*d*sin_theta);
    var_3       *= var_3;

    return var_1*var_2/var_3;
}

double Double_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a,
                                                 double d, double lambda)
{
    float sin_theta = x/sqrt(x*x + z*z);
    a = a/lambda;
    d = d/lambda;

    float var_1  = Sinc_Double(a*sin_theta);
    var_1       *= var_1;

    float var_2  = sin(TWO_PI*d*sin_theta);
    var_2       *= var_2;

    float var_3  = 2.0*sin(ONE_PI*d*sin_theta);
    var_3       *= var_3;

    return var_1*var_2/var_3;
}

long double
Double_Slit_Fraunhofer_Diffraction_Long_Double(long double x, long double z,
                                               long double a, long double d,
                                               long double lambda)
{
    float sin_theta = x/sqrtl(x*x + z*z);
    a = a/lambda;
    d = d/lambda;

    float var_1  = Sinc_Long_Double(a*sin_theta);
    var_1       *= var_1;

    float var_2  = sinl(TWO_PI*d*sin_theta);
    var_2       *= var_2;

    float var_3  = 2.0*sinl(ONE_PI*d*sin_theta);
    var_3       *= var_3;

    return var_1*var_2/var_3;
}

float Single_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a)
{
    float result = Sinc_Float(a*x/z);
    return result*result;
}

double Single_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a)
{
    double result = Sinc_Double(a*x/z);
    return result*result;
}

long double Single_Slit_Fraunhofer_Diffraction_Long_Double(long double x,
                                                           long double z,
                                                           long double a)
{
    long double result = Sinc_Long_Double(a*x/z);
    return result*result;
}
