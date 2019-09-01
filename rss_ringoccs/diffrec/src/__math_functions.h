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


extern float Max_Float(float *arr, long n_elements);

extern double Max_Double(double *arr, long n_elements);

extern long double Max_Long_Double(long double *arr, long n_elements);

extern short Max_Short(short *arr, long n_elements);

extern int Max_Int(int *arr, long n_elements);

extern long Max_Long(long *arr, long n_elements);

extern long long Max_Long_Long(long long *arr, long n_elements);


extern float Min_Float(float *arr, long n_elements);

extern double Min_Double(double *arr, long n_elements);

extern long double Min_Long_Double(long double *arr, long n_elements);

extern short Min_Short(short *arr, long n_elements);

extern int Min_Int(int *arr, long n_elements);

extern long Min_Long(long *arr, long n_elements);

extern long long Min_Long_Long(long long *arr, long n_elements);


extern float Normeq_Float(float *w_func, long n_elements);

extern double Normeq_Double(double *w_func, long n_elements);

extern long double Normeq_Long_Double(long double *w_func, long n_elements);

extern double Normeq_Short(short *w_func, long n_elements);

extern double Normeq_Int(int *w_func, long n_elements);

extern double Normeq_Long(long *w_func, long n_elements);

extern double Normeq_Long_Long(long long *w_func, long n_elements);


extern float Sinc_Float(float x);

extern double Sinc_Double(double x);

extern long double Sinc_Long_Double(long double x);

#endif