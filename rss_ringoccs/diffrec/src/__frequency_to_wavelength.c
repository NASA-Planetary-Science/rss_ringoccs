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