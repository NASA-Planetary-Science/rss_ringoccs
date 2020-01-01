#ifndef RSS_RINGOCCS_WINDOW_FUNCTIONS_H
#define RSS_RINGOCCS_WINDOW_FUNCTIONS_H

#include <complex.h>
#include "__math_functions.h"

/*  Kaiser-Bessel function with alpha = 2pi                                   */
extern float        Kaiser_Bessel_2_0_Float(float x, float W);
extern double       Kaiser_Bessel_2_0_Double(double x, double W);
extern long double  Kaiser_Bessel_2_0_Long_Double(long double x, long double W);
extern double       Kaiser_Bessel_2_0_Char(char x, double W);
extern double       Kaiser_Bessel_2_0_UChar(unsigned char x, double W);
extern double       Kaiser_Bessel_2_0_Short(short x, double W);
extern double       Kaiser_Bessel_2_0_UShort(unsigned short x, double W);
extern double       Kaiser_Bessel_2_0_Int(int x, double W);
extern double       Kaiser_Bessel_2_0_UInt(unsigned int x, double W);
extern double       Kaiser_Bessel_2_0_Long(long x, double W);
extern double       Kaiser_Bessel_2_0_ULong(unsigned long x, double W);
extern double       Kaiser_Bessel_2_0_Long_Long(long long x, double W);

extern double
Kaiser_Bessel_2_0_ULong_Long(unsigned long long x, double W);

/*  Kaiser-Bessel function with alpha = 2.5pi                                   */
extern float        Kaiser_Bessel_2_5_Float(float x, float W);
extern double       Kaiser_Bessel_2_5_Double(double x, double W);
extern long double  Kaiser_Bessel_2_5_Long_Double(long double x, long double W);
extern double       Kaiser_Bessel_2_5_Char(char x, double W);
extern double       Kaiser_Bessel_2_5_UChar(unsigned char x, double W);
extern double       Kaiser_Bessel_2_5_Short(short x, double W);
extern double       Kaiser_Bessel_2_5_UShort(unsigned short x, double W);
extern double       Kaiser_Bessel_2_5_Int(int x, double W);
extern double       Kaiser_Bessel_2_5_UInt(unsigned int x, double W);
extern double       Kaiser_Bessel_2_5_Long(long x, double W);
extern double       Kaiser_Bessel_2_5_ULong(unsigned long x, double W);
extern double       Kaiser_Bessel_2_5_Long_Long(long long x, double W);

extern double
Kaiser_Bessel_2_5_ULong_Long(unsigned long long x, double W);

/*  Kaiser-Bessel function with alpha = 3.5pi.                                */
extern float        Kaiser_Bessel_3_5_Float(float x, float W);
extern double       Kaiser_Bessel_3_5_Double(double x, double W);
extern long double  Kaiser_Bessel_3_5_Long_Double(long double x, long double W);
extern double       Kaiser_Bessel_3_5_Char(char x, double W);
extern double       Kaiser_Bessel_3_5_UChar(unsigned char x, double W);
extern double       Kaiser_Bessel_3_5_Short(short x, double W);
extern double       Kaiser_Bessel_3_5_UShort(unsigned short x, double W);
extern double       Kaiser_Bessel_3_5_Int(int x, double W);
extern double       Kaiser_Bessel_3_5_UInt(unsigned int x, double W);
extern double       Kaiser_Bessel_3_5_Long(long x, double W);
extern double       Kaiser_Bessel_3_5_ULong(unsigned long x, double W);
extern double       Kaiser_Bessel_3_5_Long_Long(long long x, double W);

extern double
Kaiser_Bessel_3_5_ULong_Long(unsigned long long x, double W);

/*  Modified Kaiser-Bessel with alpha = 2.0pi                                 */
extern float  Modified_Kaiser_Bessel_2_0_Float(float x, float W);
extern double Modified_Kaiser_Bessel_2_0_Double(double x, double W);

extern long double 
Modified_Kaiser_Bessel_2_0_Long_Double(long double x, long double W);

extern double Modified_Kaiser_Bessel_2_0_Char(char x, double W);
extern double Modified_Kaiser_Bessel_2_0_UChar(unsigned char x, double W);
extern double Modified_Kaiser_Bessel_2_0_Short(short x, double W);
extern double Modified_Kaiser_Bessel_2_0_UShort(unsigned short x, double W);
extern double Modified_Kaiser_Bessel_2_0_Int(int x, double W);
extern double Modified_Kaiser_Bessel_2_0_UInt(unsigned int x, double W);
extern double Modified_Kaiser_Bessel_2_0_Long(long x, double W);
extern double Modified_Kaiser_Bessel_2_0_ULong(unsigned long x, double W);
extern double Modified_Kaiser_Bessel_2_0_Long_Long(long long x, double W);

extern double
Modified_Kaiser_Bessel_2_0_ULong_Long(unsigned long long x, double W);

/*  Modified Kaiser-Bessel with alpha = 2.5pi                                 */
extern float  Modified_Kaiser_Bessel_2_5_Float(float x, float W);
extern double Modified_Kaiser_Bessel_2_5_Double(double x, double W);

extern long double 
Modified_Kaiser_Bessel_2_5_Long_Double(long double x, long double W);

extern double Modified_Kaiser_Bessel_2_5_Char(char x, double W);
extern double Modified_Kaiser_Bessel_2_5_UChar(unsigned char x, double W);
extern double Modified_Kaiser_Bessel_2_5_Short(short x, double W);
extern double Modified_Kaiser_Bessel_2_5_UShort(unsigned short x, double W);
extern double Modified_Kaiser_Bessel_2_5_Int(int x, double W);
extern double Modified_Kaiser_Bessel_2_5_UInt(unsigned int x, double W);
extern double Modified_Kaiser_Bessel_2_5_Long(long x, double W);
extern double Modified_Kaiser_Bessel_2_5_ULong(unsigned long x, double W);
extern double Modified_Kaiser_Bessel_2_5_Long_Long(long long x, double W);

extern double
Modified_Kaiser_Bessel_2_5_ULong_Long(unsigned long long x, double W);

/*  Modified Kaiser-Bessel with alpha = 3.5pi                                 */
extern float  Modified_Kaiser_Bessel_3_5_Float(float x, float W);
extern double Modified_Kaiser_Bessel_3_5_Double(double x, double W);

extern long double 
Modified_Kaiser_Bessel_3_5_Long_Double(long double x, long double W);

extern double Modified_Kaiser_Bessel_3_5_Char(char x, double W);
extern double Modified_Kaiser_Bessel_3_5_UChar(unsigned char x, double W);
extern double Modified_Kaiser_Bessel_3_5_Short(short x, double W);
extern double Modified_Kaiser_Bessel_3_5_UShort(unsigned short x, double W);
extern double Modified_Kaiser_Bessel_3_5_Int(int x, double W);
extern double Modified_Kaiser_Bessel_3_5_UInt(unsigned int x, double W);
extern double Modified_Kaiser_Bessel_3_5_Long(long x, double W);
extern double Modified_Kaiser_Bessel_3_5_ULong(unsigned long x, double W);
extern double Modified_Kaiser_Bessel_3_5_Long_Long(long long x, double W);

extern double
Modified_Kaiser_Bessel_3_5_ULong_Long(unsigned long long x, double W);

/*  Rectangular Window Function                                               */
extern float        Rect_Window_Float(float x, float W);
extern double       Rect_Window_Double(double x, double W);
extern long double  Rect_Window_Long_Double(long double x, long double W);
extern double       Rect_Window_Char(char x, double W);
extern double       Rect_Window_UChar(unsigned char x, double W);
extern double       Rect_Window_Short(short x, double W);
extern double       Rect_Window_UShort(unsigned short x, double W);
extern double       Rect_Window_Int(int x, double W);
extern double       Rect_Window_UInt(unsigned int x, double W);
extern double       Rect_Window_Long(long x, double W);
extern double       Rect_Window_ULong(unsigned long x, double W);
extern double       Rect_Window_Long_Long(long long x, double W);
extern double       Rect_Window_ULong_Long(unsigned long long x, double W);

/*  Squared Cosine Window Function                                            */
extern float        Coss_Window_Float(float x, float W);
extern double       Coss_Window_Double(double x, double W);
extern long double  Coss_Window_Long_Double(long double x, long double W);
extern double       Coss_Window_Char(char x, double W);
extern double       Coss_Window_UChar(unsigned char x, double W);
extern double       Coss_Window_Short(short x, double W);
extern double       Coss_Window_UShort(unsigned short x, double W);
extern double       Coss_Window_Int(int x, double W);
extern double       Coss_Window_UInt(unsigned int x, double W);
extern double       Coss_Window_Long(long x, double W);
extern double       Coss_Window_ULong(unsigned long x, double W);
extern double       Coss_Window_Long_Long(long long x, double W);
extern double       Coss_Window_ULong_Long(unsigned long long x, double W);

/*  Window Normalization Functions                                            */
extern float
Window_Normalization_Float(float *ker, long dim, float dx, float f_scale);

extern double
Window_Normalization_Double(double *ker, long dim, double dx, double f_scale);

extern long double
Window_Normalization_Long_Double(long double *ker, long dim,
                                 long double dx, long double f_scale);

extern double
Window_Normalization_Short(short *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Int(int *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Long(long *ker, long dim, double dx, double f_scale);

extern double
Window_Normalization_Long_Long(long long *ker, long dim,
                               double dx, double f_scale);

extern float
Window_Normalization_Complex_Float(complex float *ker, long dim,
                                   float dx, float f_scale);

extern double
Window_Normalization_Complex_Double(complex double *ker, long dim,
                                    double dx, double f_scale);

extern long double
Window_Normalization_Complex_Long_Double(complex long double *ker, long dim,
                                         long double dx, long double f_scale);

/*  Kaiser-Bessel function with arbitrary alpha.                              */
extern float    Kaiser_Bessel_Al_Float(float x, float W, float alpha);
extern double   Kaiser_Bessel_Al_Double(double x, double W, double alpha);

extern long double
Kaiser_Bessel_Al_Long_Double(long double x, long double W, long double alpha);

extern double Kaiser_Bessel_Al_Char(char x, double W, double alpha);
extern double Kaiser_Bessel_Al_UChar(unsigned char x, double W, double alpha);
extern double Kaiser_Bessel_Al_Short(short x, double W, double alpha);
extern double Kaiser_Bessel_Al_UShort(unsigned short x, double W, double alpha);
extern double Kaiser_Bessel_Al_Int(int x, double W, double alpha);
extern double Kaiser_Bessel_Al_UInt(unsigned int x, double W, double alpha);
extern double Kaiser_Bessel_Al_Long(long x, double W, double alpha);
extern double Kaiser_Bessel_Al_ULong(unsigned long x, double W, double alpha);
extern double Kaiser_Bessel_Al_Long_Long(long long x, double W, double alpha);

extern double
Kaiser_Bessel_Al_ULong_Long(unsigned long long x, double W, double alpha);

/*  Modified Kaiser-Bessel function with arbitrary alpha.                     */
extern float
Modified_Kaiser_Bessel_Al_Float(float x, float W, float alpha);

extern double
Modified_Kaiser_Bessel_Al_Double(double x, double W, double alpha);

extern long double
Modified_Kaiser_Bessel_Al_Long_Double(long double x, long double W,
                                      long double alpha);

extern double
Modified_Kaiser_Bessel_Al_Char(char x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_UChar(unsigned char x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_Short(short x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_UShort(unsigned short x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_Int(int x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_UInt(unsigned int x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_Long(long x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_ULong(unsigned long x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_Long_Long(long long x, double W, double alpha);

extern double
Modified_Kaiser_Bessel_Al_ULong_Long(unsigned long long x,
                                     double W, double alpha);

#endif