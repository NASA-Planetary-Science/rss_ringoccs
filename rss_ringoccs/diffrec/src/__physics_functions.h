#ifndef RSS_RINGOCCS_PHYSICS_FUNCTIONS_H
#define RSS_RINGOCCS_PHYSICS_FUNCTIONS_H

#include <math.h>
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

/*  The Fresnel Scale.                                                        */

float Fresnel_Scale_Float(float lambda, float d, float phi, float b){
    float cb_2_sp_2, sb_2;

    cb_2_sp_2  = cosf(b)*sinf(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = sinf(b);
    sb_2 *= sb_2;

    return sqrt(0.5 * lambda * d * (1.0 - cb_2_sp_2) / sb_2);
}

double Fresnel_Scale_Double(double lambda, double d, double phi, double b){
    double cb_2_sp_2, sb_2;

    cb_2_sp_2  = cos(b)*sin(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = sin(b);
    sb_2 *= sb_2;

    return sqrt(0.5 * lambda * d * (1.0 - cb_2_sp_2) / sb_2);
}


long double Fresnel_Scale_Long_Double(long double lambda, long double d,
                                      long double phi, long double b){
    long double cb_2_sp_2, sb_2;

    cb_2_sp_2  = cosl(b)*sinl(phi);
    cb_2_sp_2 *= cb_2_sp_2;

    sb_2  = sinl(b);
    sb_2 *= sb_2;

    return sqrt(0.5 * lambda * d * (1.0 - cb_2_sp_2) / sb_2);
}



#endif