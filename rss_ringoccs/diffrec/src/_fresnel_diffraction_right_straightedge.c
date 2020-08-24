#include "special_functions.h"
#include "_math_constants.h"

/**********Right Straightedge Diffraction Using Fresnel Approximation**********/

complex float Right_Straightedge_Diffraction_Float(float x, float edge,
                                                   float F){
    complex float T_hat;
    float re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sin_Float(x);
    re = Fresnel_Cos_Float(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5-SQRT_2_BY_PI*T_hat;
}

complex double Right_Straightedge_Diffraction_Double(double x, double edge,
                                                     double F){
    complex double T_hat;
    double re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sin_Double(x);
    re = Fresnel_Cos_Double(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5-SQRT_2_BY_PI*T_hat;
}

complex long double Right_Straightedge_Diffraction_Long_Double(long double x,
                                                               long double edge,
                                                               long double F){
    complex long double T_hat;
    long double re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sin_Long_Double(x);
    re = Fresnel_Cos_Long_Double(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5-SQRT_2_BY_PI*T_hat;
}

/*  For all integer types, convert to double and compute.                     */
double Right_Straightedge_Diffraction_Char(char x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_UChar(unsigned char x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_Short(short x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_UShort(unsigned short x, double a,
                                             double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_Int(int x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_UInt(unsigned int x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_Long(long x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_ULong(unsigned long x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_Long_Long(long long x, double a, double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}

double Right_Straightedge_Diffraction_ULong_Long(unsigned long long x, double a,
                                                 double F)
{
    return Right_Straightedge_Diffraction_Double((double)x, a, F);
}
