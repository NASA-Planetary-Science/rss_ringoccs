#include <rss_ringoccs_special_functions.h>
#include <rss_ringoccs_math_constants.h>

/***********Left Straightedge Diffraction Using Fresnel Approximation**********/

complex float Left_Straightedge_Diffraction_Float(float x, float edge, float F)
{
    complex float T_hat;
    float re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sin_Float(x);
    re = Fresnel_Cos_Float(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5 + SQRT_2_BY_PI*T_hat;
}

complex double Left_Straightedge_Diffraction_Double(double x, double edge,
                                                    double F)
{
    complex double T_hat;
    double re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sin_Double(x);
    re = Fresnel_Cos_Double(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5 + SQRT_2_BY_PI*T_hat;
}

complex long double Left_Straightedge_Diffraction_Long_Double(
    long double x, long double edge, long double F)
{
    complex long double T_hat;
    long double re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sin_Long_Double(x);
    re = Fresnel_Cos_Long_Double(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5 + SQRT_2_BY_PI*T_hat;
}
