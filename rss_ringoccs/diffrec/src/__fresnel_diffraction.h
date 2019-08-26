#ifndef RSS_RINGOCCS_FRESNEL_DIFFRACTION_H
#define RSS_RINGOCCS_FRESNEL_DIFFRACTION_H

#include "__fresnel_integrals.h"
#include "__math_constants.h"

/*******************************************************************************
 *------------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.   *
 * They are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.                *
 ******************************************************************************/

/*************Square Well Diffraction Using Fresnel Approximation**************/

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Solution_Float                                *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from a plane wave incident on a       *
 *      square well, assuming the Fresnel approximation is valid.             *
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
complex float Square_Well_Diffraction_Solution_Float(float x, float a,
                                                     float b, float F){
    float arg1 = SQRT_PI_BY_2*(a-x)/F;
    float arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex float result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
                           Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return 1.0 - (0.5 - 0.5*_Complex_I)*result;
}

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Solution_Double                               *
 *  Purpose:                                                                  *
 *      Same as Square_Well_Diffraction_Solution_Float, but for doubles.      *
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
complex double Square_Well_Diffraction_Solution_Double(double x, double a,
                                                       double b, double F){
    double arg1 = SQRT_PI_BY_2*(a-x)/F;
    double arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex double result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
                            Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return 1.0 - (0.5 - 0.5*_Complex_I)*result;
}

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Solution_Long_Double                          *
 *  Purpose:                                                                  *
 *      Same as Square_Well_Diffraction_Solution_Float, but for long doubles. *
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
complex long double Square_Well_Diffraction_Solution_Long_Double(double x,
                                                                 double a,
                                                                 double b,
                                                                 double F){
    long double arg1 = SQRT_PI_BY_2*(a-x)/F;
    long double arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex long double result;
    
    result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
             Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return 1.0 - (0.5 - 0.5*_Complex_I)*result;
}

/*--------Inverted Square Well Diffraction Using Fresnel Approximation--------*/

/******************************************************************************
 *  Function:                                                                 *
 *      Inverted_Square_Well_Diffraction_Solution_Float                       *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from a plane wave incident on a       *
 *      square well, assuming the Fresnel approximation is valid.             *
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
complex float Inverted_Square_Well_Diffraction_Solution_Float(float x, float a,
                                                              float b, float F){
    float arg1 = SQRT_PI_BY_2*(a-x)/F;
    float arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex float result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
                           Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return (0.5 - 0.5*_Complex_I)*result;
}

/******************************************************************************
 *  Function:                                                                 *
 *      Inverted_Square_Well_Diffraction_Solution_Double                      *
 *  Purpose:                                                                  *
 *      Acts as Inverted_Square_Well_Diffraction_Solution_Float for doubles.  *
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
complex double Inverted_Square_Well_Diffraction_Solution_Double(double x,
                                                                double a,
                                                                double b,
                                                                double F){
    double arg1 = SQRT_PI_BY_2*(a-x)/F;
    double arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex double result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
                            Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return (0.5 - 0.5*_Complex_I)*result;
}

/******************************************************************************
 *  Function:                                                                 *
 *      Inverted_Square_Well_Diffraction_Solution_Long_Double                 *
 *  Purpose:                                                                  *
 *      Inverted_Square_Well_Diffraction_Solution_Float for long doubles.     *
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
complex long double Inverted_Square_Well_Diffraction_Solution_Long_Double(
    double x, double a, double b, double F
){
    long double arg1 = SQRT_PI_BY_2*(a-x)/F;
    long double arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex long double result;
    result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
             Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return (0.5 - 0.5*_Complex_I)*result;
}

/*--------Inverted Square Well Diffraction Using Fresnel Approximation--------*/

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
        (Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(a)   - 
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(b) +
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(b) -
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(a) +
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(a));
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
        (Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(a)   - 
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(b) +
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(b) -
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(a) +
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(a));
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
        (Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(a)   - 
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(b) +
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(a));

    re = 1.0 - SQRT_ONE_BY_2_PI *
        (Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(b) -
         Fresnel_Cosine_Heald_Rational_EPS_Minus_Eight_Func(a) +
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(b)   -
         Fresnel_Sine_Heald_Rational_EPS_Minus_Eight_Func(a));
    return atan2(im, re);
}

#endif