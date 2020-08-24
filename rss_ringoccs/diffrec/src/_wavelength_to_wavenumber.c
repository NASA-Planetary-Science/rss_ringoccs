/*  Coefficients for the Taylor series and asymptotic expansions found here.  */
#include "_math_constants.h"

/*  Main header for the math functions. Contains <math.h> as well.            */
#include "special_functions.h"

float Wavelength_To_Wavenumber_Float(float lambda){
    return TWO_PI / lambda;
}

double Wavelength_To_Wavenumber_Double(double lambda){
    return TWO_PI / lambda;
}

long double Wavelength_To_Wavenumber_Long_Double(long double lambda){
    return TWO_PI / lambda;
}

/*  For all integer types, convert to double and compute.                     */
RSSRINGOCCSNonFloatInputForFloatOutput(Wavelength_To_Wavenumber);
