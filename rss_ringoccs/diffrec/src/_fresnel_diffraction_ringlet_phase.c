#include "special_functions.h"
#include "_math_constants.h"

/***************Ringlet Diffraction Using Fresnel Approximation****************/

/******************************************************************************
 *  Function:                                                                 *
 *      Ringlet_Diffraction_Phase_Float                                       *
 *  Purpose:                                                                  *
 *      Compute the phase from the diffraction pattern of a square well.      *
 *  Arguments:                                                                *
 *      x (float):                                                            *
 *          The location on the x-axis for the point being computed.          *
 *      a (float):                                                            *
 *          The left-most endpoint of the ringlet.                            *
 *      b (float):                                                            *
 *          The right-most endpoint of the ringlet.                           *
 *      F (float):                                                            *
 *          The Fresnel scale.                                                *
 *  Notes:                                                                    *
 *      1.) This function relies on the C99 standard, or higher.              *
 ******************************************************************************/
float Ringlet_Diffraction_Phase_Float(float x, float a, float b, float F)
{
    float re, im;

    a = SQRT_PI_BY_2*(a-x)/F;
    b = SQRT_PI_BY_2*(b-x)/F;

    im = SQRT_ONE_BY_2_PI *
        (Fresnel_Sin_Float(b) - Fresnel_Sin_Float(a) - 
         Fresnel_Cos_Float(b) + Fresnel_Cos_Float(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cos_Float(b) - Fresnel_Cos_Float(a) +
         Fresnel_Sin_Float(b) - Fresnel_Sin_Float(a));
    return atan2(im, re);
}

/******************************************************************************
 *  Function:                                                                 *
 *      Ringlet_Diffraction_Phase_Double                                      *
 *  Purpose:                                                                  *
 *      Compute the phase from the diffraction pattern of a square well.      *
 *  Arguments:                                                                *
 *      x (float):                                                            *
 *          The location on the x-axis for the point being computed.          *
 *      a (float):                                                            *
 *          The left-most endpoint of the ringlet.                            *
 *      b (float):                                                            *
 *          The right-most endpoint of the ringlet.                           *
 *      F (float):                                                            *
 *          The Fresnel scale.                                                *
 *  Notes:                                                                    *
 *      1.) This function relies on the C99 standard, or higher.              *
 ******************************************************************************/
double Ringlet_Diffraction_Phase_Double(double x, double a, double b, double F)
{
    double re, im;

    a = SQRT_PI_BY_2*(a-x)/F;
    b = SQRT_PI_BY_2*(b-x)/F;

    im = SQRT_ONE_BY_2_PI *
        (Fresnel_Sin_Double(b) - Fresnel_Sin_Double(a) - 
         Fresnel_Cos_Double(b) + Fresnel_Cos_Double(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cos_Double(b) - Fresnel_Cos_Double(a) +
         Fresnel_Sin_Double(b) - Fresnel_Sin_Double(a));
    return atan2(im, re);
}

/******************************************************************************
 *  Function:                                                                 *
 *      Ringlet_Diffraction_Phase_Long_Double                                 *
 *  Purpose:                                                                  *
 *      Compute the phase from the diffraction pattern of a square well.      *
 *  Arguments:                                                                *
 *      x (float):                                                            *
 *          The location on the x-axis for the point being computed.          *
 *      a (float):                                                            *
 *          The left-most endpoint of the ringlet.                            *
 *      b (float):                                                            *
 *          The right-most endpoint of the ringlet.                           *
 *      F (float):                                                            *
 *          The Fresnel scale.                                                *
 *  Notes:                                                                    *
 *      1.) This function relies on the C99 standard, or higher.              *
 ******************************************************************************/
long double Ringlet_Diffraction_Phase_Long_Double(long double x, long double a,
                                                  long double b, long double F)
{
    long double re, im;

    a = SQRT_PI_BY_2*(a-x)/F;
    b = SQRT_PI_BY_2*(b-x)/F;

    im = SQRT_ONE_BY_2_PI *
        (Fresnel_Sin_Long_Double(b)   -
         Fresnel_Sin_Long_Double(a)   - 
         Fresnel_Cos_Long_Double(b) +
         Fresnel_Cos_Long_Double(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cos_Long_Double(b) - Fresnel_Cos_Long_Double(a) +
         Fresnel_Sin_Long_Double(b) - Fresnel_Sin_Long_Double(a));
    return atan2(im, re);
}
