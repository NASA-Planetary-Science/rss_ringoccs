#ifndef RSS_RINGOCCS_PHYSICS_FUNCTIONS_H
#define RSS_RINGOCCS_PHYSICS_FUNCTIONS_H

#include "__math_constants.h"

/*  Wavelength to Wavenumber function.                                        */

/* Convert wavelength to wavenumber for a float input. Returns float.         */
float Wavelength_To_Wavenumber_Float(float lambda){
    return TWO_PI / lambda;
}

double Wavelength_To_Wavenumber_Double(double lambda){
    return TWO_PI / lambda;
}

long double Wavelength_To_Wavenumber_Long_Double(long double lambda){
    return TWO_PI / lambda;
}

/*  Frequency to Wavelength function.                                         */
float Frequency_To_Wavelength_Float(float frequency){
    return SPEED_OF_LIGHT_KM / frequency;
}

double Frequency_To_Wavelength_Double(double frequency){
    return SPEED_OF_LIGHT_KM / frequency;
}

long double Frequency_To_Wavelength_Long_Double(long double frequency){
    return SPEED_OF_LIGHT_KM / frequency;
}

#endif