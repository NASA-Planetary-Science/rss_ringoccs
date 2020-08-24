/*  Coefficients for the Taylor series and asymptotic expansions found here.  */
#include "_math_constants.h"

/*  Main header for the math functions. Contains <math.h> as well.            */
#include "special_functions.h"

float Frequency_To_Wavelength_Float(float frequency){
    return SPEED_OF_LIGHT_KM / frequency;
}

double Frequency_To_Wavelength_Double(double frequency){
    return SPEED_OF_LIGHT_KM / frequency;
}

long double Frequency_To_Wavelength_Long_Double(long double frequency){
    return SPEED_OF_LIGHT_KM / frequency;
}

/*  For all integer types, convert to double and compute.                     */
RSSRINGOCCSNonFloatInputForFloatOutput(Frequency_To_Wavelength);
