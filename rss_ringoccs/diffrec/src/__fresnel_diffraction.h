#ifndef RSS_RINGOCCS_FRESNEL_DIFFRACTION_H
#define RSS_RINGOCCS_FRESNEL_DIFFRACTION_H

#include "__fresnel_integrals.h"

/*******************************************************************************
 *------------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.   *
 * They are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.                *
 ******************************************************************************/

/*************Square Well Diffraction Using Fresnel Approximation**************/

extern complex float Square_Well_Diffraction_Float(float x, float a,
                                                   float b, float F);

extern complex double Square_Well_Diffraction_Double(double x, double a,
                                                     double b, double F);

extern complex long double Square_Well_Diffraction_Long_Double(long double x,
                                                               long double a,
                                                               long double b,
                                                               long double F);

/*********nverted Square Well Diffraction Using Fresnel Approximation**********/

extern complex float Inverted_Square_Well_Diffraction_Float(float x, float a,
                                                            float b, float F);

extern complex double Inverted_Square_Well_Diffraction_Double(double x,
                                                              double a,
                                                              double b,
                                                              double F);

extern complex long double Inverted_Square_Well_Diffraction_Long_Double(
    long double x, long double a, long double b, long double F
);

/*--------Inverted Square Well Diffraction Using Fresnel Approximation--------*/

extern float Square_Well_Diffraction_Phase_Float(float x, float a,
                                                 float b, float F);

extern double Square_Well_Diffraction_Phase_Double(double x, double a,
                                                   double b, double F);

extern long double Square_Well_Diffraction_Phase_Long_Double(long double x,
                                                             long double a,
                                                             long double b,
                                                             long double F);

/*--------Right Straight-Edge Diffraction Using Fresnel Approximation---------*/

extern complex float Right_Straightedge_Diffraction_Float(float x, float edge,
                                                          float F);

extern complex double Right_Straightedge_Diffraction_Double(double x,
                                                            double edge,
                                                            double F);

extern complex long double Right_Straightedge_Diffraction_Long_Double(
    long double x, long double edge, long double F
);

/*---------Left Straight-Edge Diffraction Using Fresnel Approximation---------*/

extern complex float Left_Straightedge_Diffraction_Float(float x, float edge,
                                                         float F);

extern complex double Left_Straightedge_Diffraction_Double(double x,
                                                           double edge,
                                                           double F);

extern complex long double Left_Straightedge_Diffraction_Long_Double(
    long double x, long double edge, long double F
);

#endif