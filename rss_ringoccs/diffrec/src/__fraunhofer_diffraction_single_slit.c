#include "__fraunhofer_diffraction.h"

float Single_Slit_Fraunhofer_Diffraction_Float(float x, float z, float a)
{
    float result = Sinc_Float(a*x/z);
    return result*result;
}

double Single_Slit_Fraunhofer_Diffraction_Double(double x, double z, double a)
{
    double result = Sinc_Double(a*x/z);
    return result*result;
}

long double Single_Slit_Fraunhofer_Diffraction_Long_Double(long double x,
                                                           long double z,
                                                           long double a)
{
    long double result = Sinc_Long_Double(a*x/z);
    return result*result;
}