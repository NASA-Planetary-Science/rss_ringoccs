#include "__fresnel_diffraction.h"

/*************Square Well Diffraction Using Fresnel Approximation**************/

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Float                                         *
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
complex float Square_Well_Diffraction_Float(float x, float a, float b, float F){
    float arg1 = SQRT_PI_BY_2*(a-x)/F;
    float arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex float result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
                           Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return 1.0 - (0.5 - 0.5*_Complex_I)*result;
}

/******************************************************************************
 *  Function:                                                                 *
 *      Square_Well_Diffraction_Double                                        *
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
complex double Square_Well_Diffraction_Double(double x, double a,
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
 *      Square_Well_Diffraction_Long_Double                                   *
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
complex long double Square_Well_Diffraction_Long_Double(
    long double x, long double a, long double b, long double F
){

    long double arg1 = SQRT_PI_BY_2*(a-x)/F;
    long double arg2 = SQRT_PI_BY_2*(b-x)/F;
    complex long double result;
    
    result = Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg2) -
             Fresnel_Heald_Rational_EPS_Minus_Eight_Func(arg1);
    result *= SQRT_2_BY_PI;

    return 1.0 - (0.5 - 0.5*_Complex_I)*result;
}