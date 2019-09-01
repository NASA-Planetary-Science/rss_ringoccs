#ifndef RSS_RINGOCCS_WINDOW_FUNCTIONS_H
#define RSS_RINGOCCS_WINDOW_FUNCTIONS_H

#include <math.h>
#include "__math_constants.h"

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

#endif