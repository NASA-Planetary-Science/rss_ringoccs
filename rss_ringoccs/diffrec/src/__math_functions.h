#ifndef RSS_RINGOCCS_MATH_FUNCTIONS_H
#define RSS_RINGOCCS_MATH_FUNCTIONS_H

#include <math.h>
#include "__math_constants.h"

extern float LambertW_Float(float x);

extern double LambertW_Double(double x);

extern long double LambertW_Long_Double(long double x);

extern void Legendre_Polynomials(double *legendre_p, double x, int order);

extern void Alt_Legendre_Polynomials(double *poly, double *legendre_p,
                                     int order);

extern void Fresnel_Kernel_Coefficients(double *fresnel_ker_coeffs,
                                        double *legendre_p,
                                        double *alt_legendre_p,
                                        double Legendre_Coeff, int order);

extern double Resolution_Inverse_Float(float x);

extern double Resolution_Inverse_Double(double x);

extern long double Resolution_Inverse_Long_Double(long double x);

extern float BesselJ0_Float(float x);

extern double BesselJ0_Double(double x);

extern long double BesselJ0_Long_Double(long double x);


#endif