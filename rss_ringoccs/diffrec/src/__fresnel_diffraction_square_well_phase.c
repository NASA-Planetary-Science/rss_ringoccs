#include "__fresnel_diffraction.h"

/*********Inverted Square Well Diffraction Using Fresnel Approximation*********/

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Phase_Float                                   *
 *  Purpose:                                                                  *
 *      Compute the phase from the diffraction pattern of a square well.      *
 *  Arguments:                                                                *
 *      x (float):                                                            *
 *          The location on the x-axis for the point being computed.          *
 *      a (float):                                                            *
 *          The left-most endpoint of the square well.                        *
 *      b (float):                                                            *
 *          The right-most endpoint of the square well.                       *
 *      F (float):                                                            *
 *          The Fresnel scale.                                                *
 *  Notes:                                                                    *
 *      1.) This function relies on the C99 standard, or higher.              *
 ******************************************************************************/
float Square_Well_Diffraction_Phase_Float(float x, float a, float b, float F){
    float re, im;

    a = SQRT_PI_BY_2*(a-x)/F;
    b = SQRT_PI_BY_2*(b-x)/F;

    im = SQRT_ONE_BY_2_PI *
        (Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Float(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Float(a)   - 
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Float(b) +
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Float(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Float(b) -
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Float(a) +
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Float(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Float(a));
    return atan2(im, re);
}

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Phase_Double                                  *
 *  Purpose:                                                                  *
 *      Compute the phase from the diffraction pattern of a square well.      *
 *  Arguments:                                                                *
 *      x (double):                                                           *
 *          The location on the x-axis for the point being computed.          *
 *      a (double):                                                           *
 *          The left-most endpoint of the square well.                        *
 *      b (double):                                                           *
 *          The right-most endpoint of the square well.                       *
 *      F (double):                                                           *
 *          The Fresnel scale.                                                *
 *  Notes:                                                                    *
 *      1.) This function relies on the C99 standard, or higher.              *
 ******************************************************************************/
double Square_Well_Diffraction_Phase_Double(double x, double a,
                                            double b, double F){
    double re, im;

    a = SQRT_PI_BY_2*(a-x)/F;
    b = SQRT_PI_BY_2*(b-x)/F;

    im = SQRT_ONE_BY_2_PI *
        (Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Double(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Double(a)   - 
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Double(b) +
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Double(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Double(b) -
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Double(a) +
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Double(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Double(a));
    return atan2(im, re);
}

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Phase_Long_Double                             *
 *  Purpose:                                                                  *
 *      Compute the phase from the diffraction pattern of a square well.      *
 *  Arguments:                                                                *
 *      x (double):                                                           *
 *          The location on the x-axis for the point being computed.          *
 *      a (double):                                                           *
 *          The left-most endpoint of the square well.                        *
 *      b (double):                                                           *
 *          The right-most endpoint of the square well.                       *
 *      F (double):                                                           *
 *          The Fresnel scale.                                                *
 *  Notes:                                                                    *
 *      1.) This function relies on the C99 standard, or higher.              *
 ******************************************************************************/
long double Square_Well_Diffraction_Phase_Long_Double(long double x,
                                                      long double a,
                                                      long double b,
                                                      long double F){
    long double re, im;

    a = SQRT_PI_BY_2*(a-x)/F;
    b = SQRT_PI_BY_2*(b-x)/F;

    im = SQRT_ONE_BY_2_PI *
        (Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Long_Double(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Long_Double(a)   - 
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Long_Double(b) +
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Long_Double(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Long_Double(b) -
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Long_Double(a) +
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Long_Double(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Long_Double(a));
    return atan2(im, re);
}