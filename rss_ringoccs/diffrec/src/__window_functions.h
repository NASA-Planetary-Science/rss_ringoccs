#ifndef RSS_RINGOCCS_WINDOW_FUNCTIONS_H
#define RSS_RINGOCCS_WINDOW_FUNCTIONS_H

#include "__math_functions.h"

/*  Kaiser-Bessel function with alpha = 2pi                                   */
extern float Kaiser_Bessel_2_0_Float(float x, float W);

extern double Kaiser_Bessel_2_0_Double(double x, double W);

extern long double Kaiser_Bessel_2_0_Long_Double(long double x, long double W);

/*  Kaiser-Bessel function with alpha = 2.5pi.                                */
extern float Kaiser_Bessel_2_5_Float(float x, float W);

extern double Kaiser_Bessel_2_5_Double(double x, double W);

extern long double Kaiser_Bessel_2_5_Long_Double(long double x, long double W);

/*  Kaiser-Bessel function with alpha = 3.5pi.                                */
extern float Kaiser_Bessel_3_5_Float(float x, float W);

extern double Kaiser_Bessel_3_5_Double(double x, double W);

extern long double Kaiser_Bessel_3_5_Long_Double(long double x, long double W);

/*  Modified Kaiser-Bessel with alpha = 2.0pi                                 */
extern float Modified_Kaiser_Bessel_2_0_Float(float x, float W);

extern double Modified_Kaiser_Bessel_2_0_Double(double x, double W);

extern long double Modified_Kaiser_Bessel_2_0_Long_Double(long double x,
                                                          long double W);

/*  Modified Kaiser-Bessel with alpha = 2.5pi                                 */
extern float Modified_Kaiser_Bessel_2_5_Float(float x, float W);

extern double Modified_Kaiser_Bessel_2_5_Double(double x, double W);

extern long double Modified_Kaiser_Bessel_2_5_Long_Double(long double x,
                                                          long double W);

/*  Modified Kaiser-Bessel with alpha = 2.5pi                                 */
extern float Modified_Kaiser_Bessel_3_5_Float(float x, float W);

extern double Modified_Kaiser_Bessel_3_5_Double(double x, double W);

extern long double Modified_Kaiser_Bessel_3_5_Long_Double(long double x,
                                                          long double W);

/*  Rectangular Window Function                                               */
extern float Rect_Window_Float(float x, float W);

extern double Rect_Window_Double(double x, double W);

extern long double Rect_Window_Long_Double(long double x, long double W);

/*  Squared Cosine Window Function                                            */
extern float Coss_Window_Float(float x, float W);

extern double Coss_Window_Double(double x, double W);

extern long double Coss_Window_Long_Double(long double x, long double W);


/*  Window Normalization Functions                                            */
extern float Window_Normalization_Float(float *ker, long dim,
                                        float dx, float f_scale);

extern double Window_Normalization_Double(double *ker, long dim,
                                          double dx, double f_scale);

extern long double Window_Normalization_Long_Double(long double *ker, long dim,
                                                    long double dx,
                                                    long double f_scale);

extern double Window_Normalization_Short(short *ker, long dim,
                                         double dx, double f_scale);

extern double Window_Normalization_Int(int *ker, long dim,
                                       double dx, double f_scale);

extern double Window_Normalization_Long(long *ker, long dim,
                                        double dx, double f_scale);

double Window_Normalization_Long_Long(long long *ker, long dim,
                                      double dx, double f_scale);

/*  Kaiser-Bessel function with arbitrary alpha.                              */
extern float Kaiser_Bessel_Al_Float(float x, float W, float alpha);

extern double Kaiser_Bessel_Al_Double(double x, double W, double alpha);

extern long double Kaiser_Bessel_Al_Long_Double(long double x, long double W,
                                                long double alpha);

/*  Modified Kaiser-Bessel function with arbitrary alpha.                     */
extern float Modified_Kaiser_Bessel_Al_Float(float x, float W, float alpha);

extern double Modified_Kaiser_Bessel_Al_Double(double x, double W,
                                               double alpha);

extern long double Modified_Kaiser_Bessel_Al_Long_Double(long double x,
                                                         long double W,
                                                         long double alpha);

#endif