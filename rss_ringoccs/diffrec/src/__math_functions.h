#ifndef RSS_RINGOCCS_MATH_FUNCTIONS_H
#define RSS_RINGOCCS_MATH_FUNCTIONS_H

#include <math.h>
#include "__math_constants.h"
#include "__get_array.h"

typedef char bool;
#define True 1
#define False 0

extern float LambertW_Float(float x);

extern double LambertW_Double(double x);

extern long double LambertW_Long_Double(long double x);


extern void Legendre_Polynomials(double *legendre_p, double x, int order);

extern void Alt_Legendre_Polynomials(double *poly, double *legendre_p,
                                     int order);

extern void
Fresnel_Kernel_Coefficients(double *fresnel_ker_coeffs, double *legendre_p,
                            double *alt_legendre_p, double Legendre_Coeff,
                            int order, bool IsEven);


extern float Resolution_Inverse_Float(float x);

extern double Resolution_Inverse_Double(double x);

extern long double Resolution_Inverse_Long_Double(long double x);


extern float        BesselJ0_Float(float x);
extern double       BesselJ0_Double(double x);
extern long double  BesselJ0_Long_Double(long double x);
extern double       BesselJ0_Char(char x);
extern double       BesselJ0_UChar(unsigned char x);
extern double       BesselJ0_Short(short x);
extern double       BesselJ0_UShort(unsigned short x);
extern double       BesselJ0_Int(int x);
extern double       BesselJ0_UInt(unsigned int x);
extern double       BesselJ0_Long(long x);
extern double       BesselJ0_ULong(unsigned long x);
extern double       BesselJ0_Long_Long(long long x);
extern double       BesselJ0_ULong_Long(unsigned long long x);

extern float        BesselI0_Float(float x);
extern double       BesselI0_Double(double x);
extern long double  BesselI0_Long_Double(long double x);
extern double       BesselI0_Char(char x);
extern double       BesselI0_UChar(unsigned char x);
extern double       BesselI0_Short(short x);
extern double       BesselI0_UShort(unsigned short x);
extern double       BesselI0_Int(int x);
extern double       BesselI0_UInt(unsigned int x);
extern double       BesselI0_Long(long x);
extern double       BesselI0_ULong(unsigned long x);
extern double       BesselI0_Long_Long(long long x);
extern double       BesselI0_ULong_Long(unsigned long long x);


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


extern float        Sinc_Float(float x);
extern double       Sinc_Double(double x);
extern long double  Sinc_Long_Double(long double x);
extern double       Sinc_Char(char x);
extern double       Sinc_UChar(unsigned char x);
extern double       Sinc_Short(short x);
extern double       Sinc_UShort(unsigned short x);
extern double       Sinc_Int(int x);
extern double       Sinc_UInt(unsigned int x);
extern double       Sinc_Long(long x);
extern double       Sinc_ULong(unsigned long x);
extern double       Sinc_Long_Long(long long x);
extern double       Sinc_ULong_Long(unsigned long long x);

#endif