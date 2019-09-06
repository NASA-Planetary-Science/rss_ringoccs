#include "__fresnel_diffraction.h"

/***********Left Straightedge Diffraction Using Fresnel Approximation**********/

complex float Left_Straightedge_Diffraction_Float(float x, float edge, float F)
{
    complex float T_hat;
    float re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sine_Taylor_to_Asymptotic_Float(x);
    re = Fresnel_Cosine_Taylor_to_Asymptotic_Float(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5 + SQRT_2_BY_PI*T_hat;
}

complex double Left_Straightedge_Diffraction_Double(double x, double edge,
                                                    double F)
{
    complex double T_hat;
    double re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sine_Taylor_to_Asymptotic_Double(x);
    re = Fresnel_Cosine_Taylor_to_Asymptotic_Double(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5 + SQRT_2_BY_PI*T_hat;
}

complex long double Left_Straightedge_Diffraction_Long_Double(
    long double x, long double edge, long double F)
{
    complex long double T_hat;
    long double re, im;

    x = SQRT_PI_BY_2*(edge-x)/F;

    im = Fresnel_Sine_Taylor_to_Asymptotic_Long_Double(x);
    re = Fresnel_Cosine_Taylor_to_Asymptotic_Long_Double(x);

    T_hat = (0.5 - 0.5*_Complex_I)*(re+_Complex_I*im);
    return 0.5 + SQRT_2_BY_PI*T_hat;
}