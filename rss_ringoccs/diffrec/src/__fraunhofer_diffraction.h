#ifndef RSS_RINGOCCS_FRAUNHOFER_DIFFRACTION_H
#define RSS_RINGOCCS_FRAUNHOFER_DIFFRACTION_H

#include <math.h>
#include "__sinc.h"
#include "__math_constants.h"

float Single_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a){
    float result = Sinc_Float(a*x/z);
    return result*result;
}

double Single_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a){
    double result = Sinc_Double(a*x/z);
    return result*result;
}

long double Single_Slit_Fraunhofer_Diffraction_Long_Double(long double x,
                                                           long double z,
                                                           long double a){
    long double result = Sinc_Long_Double(a*x/z);
    return result*result;
}

float Double_Slit_Fraunhofer_Diffraction_Float(float x, float z,
                                               float a, float d){
    float var_1, var_2, var_3;

    var_1  = Sinc_Float(a*x/z);
    var_1 *= var_1;

    var_2  = sinf(TWO_PI*d*x/z);
    var_2 *= var_2;

    var_3  = 2.0*sinf(ONE_PI*d*x/z);
    var_3 *= var_3;

    return var_1*var_2/var_3;
}

double Double_Slit_Fraunhofer_Diffraction_Double(double x, double z,
                                                 double a, double d){
    double var_1, var_2, var_3;

    var_1  = Sinc_Double(a*x/z);
    var_1 *= var_1;

    var_2  = sin(TWO_PI*d*x/z);
    var_2 *= var_2;

    var_3  = 2.0*sin(ONE_PI*d*x/z);
    var_3 *= var_3;

    return var_1*var_2/var_3;
}

long double Double_Slit_Fraunhofer_Diffraction_Long_Double(long double x,
                                                           long double z,
                                                           long double a,
                                                           long double d){
    long double var_1, var_2, var_3;

    var_1  = Sinc_Long_Double(a*x/z);
    var_1 *= var_1;

    var_2  = sinl(TWO_PI*d*x/z);
    var_2 *= var_2;

    var_3  = 2.0*sinl(ONE_PI*d*x/z);
    var_3 *= var_3;

    return var_1*var_2/var_3;
}

#endif