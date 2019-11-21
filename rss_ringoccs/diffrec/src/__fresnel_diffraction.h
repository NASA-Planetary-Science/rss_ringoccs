#ifndef RSS_RINGOCCS_FRESNEL_DIFFRACTION_H
#define RSS_RINGOCCS_FRESNEL_DIFFRACTION_H

#include "__fresnel_integrals.h"

/*******************************************************************************
 *------------------------------DEFINE C FUNCTIONS-----------------------------*
 * These are functions written in pure C without the use of the Numpy-C API.   *
 * They are used to define various special functions. They will be wrapped in  *
 * a form that is useable with the Python interpreter later on.                *
 ******************************************************************************/

/*-------------Ringlet Diffraction Using Fresnel Approximation----------------*/

extern complex float
Ringlet_Diffraction_Float(float x, float a, float b, float F);

extern complex double
Ringlet_Diffraction_Double(double x, double a, double b, double F);

extern complex long double
Ringlet_Diffraction_Long_Double(long double x, long double a,
                                long double b, long double F);

extern double
Ringlet_Diffraction_Char(char x, double a, double b, double F);

extern double
Ringlet_Diffraction_UChar(unsigned char x, double a, double b, double F);

extern double
Ringlet_Diffraction_Short(short x, double a, double b, double F);

extern double
Ringlet_Diffraction_UShort(unsigned short x, double a, double b, double F);

extern double
Ringlet_Diffraction_Int(int x, double a, double b, double F);

extern double
Ringlet_Diffraction_UInt(unsigned int x, double a, double b, double F);

extern double
Ringlet_Diffraction_Long(long x, double a, double b, double F);

extern double
Ringlet_Diffraction_ULong(unsigned long x, double a, double b, double F);

extern double
Ringlet_Diffraction_Long_Long(long long x, double a, double b, double F);

extern double
Ringlet_Diffraction_ULong_Long(unsigned long long x, double a,
                               double b, double F);

/*----------------Gap Diffraction Using Fresnel Approximation-----------------*/

extern complex float
Gap_Diffraction_Float(float x, float a, float b, float F);

extern complex double
Gap_Diffraction_Double(double x, double a, double b, double F);

extern complex long double
Gap_Diffraction_Long_Double(long double x, long double a,
                            long double b, long double F);

extern double
Gap_Diffraction_Char(char x, double a, double b, double F);

extern double
Gap_Diffraction_UChar(unsigned char x, double a, double b, double F);

extern double
Gap_Diffraction_Short(short x, double a, double b, double F);

extern double
Gap_Diffraction_UShort(unsigned short x, double a, double b, double F);

extern double
Gap_Diffraction_Int(int x, double a, double b, double F);

extern double
Gap_Diffraction_UInt(unsigned int x, double a, double b, double F);

extern double
Gap_Diffraction_Long(long x, double a, double b, double F);

extern double
Gap_Diffraction_ULong(unsigned long x, double a, double b, double F);

extern double
Gap_Diffraction_Long_Long(long long x, double a, double b, double F);

extern double
Gap_Diffraction_ULong_Long(unsigned long long x, double a,
                           double b, double F);

/*-----------Ringlet Diffraction Phase Using Fresnel Approximation------------*/

extern float
Ringlet_Diffraction_Phase_Float(float x, float a, float b, float F);

extern double
Ringlet_Diffraction_Phase_Double(double x, double a, double b, double F);

extern long double
Ringlet_Diffraction_Phase_Long_Double(long double x, long double a,
                                      long double b, long double F);

/*--------Right Straight-Edge Diffraction Using Fresnel Approximation---------*/

extern complex float
Right_Straightedge_Diffraction_Float(float x, float edge, float F);

extern complex double
Right_Straightedge_Diffraction_Double(double x, double edge, double F);

extern complex long double
Right_Straightedge_Diffraction_Long_Double(long double x, long double edge,
                                           long double F);

extern double
Right_Straightedge_Diffraction_Char(char x, double a, double F);

extern double
Right_Straightedge_Diffraction_UChar(unsigned char x, double a, double F);

extern double
Right_Straightedge_Diffraction_Short(short x, double a, double F);

extern double
Right_Straightedge_Diffraction_UShort(unsigned short x, double a, double F);

extern double
Right_Straightedge_Diffraction_Int(int x, double a, double F);

extern double
Right_Straightedge_Diffraction_UInt(unsigned int x, double a, double F);

extern double
Right_Straightedge_Diffraction_Long(long x, double a, double F);

extern double
Right_Straightedge_Diffraction_ULong(unsigned long x, double a, double F);

extern double
Right_Straightedge_Diffraction_Long_Long(long long x, double a, double F);

extern double
Right_Straightedge_Diffraction_ULong_Long(unsigned long long x,
                                          double a, double F);

/*---------Left Straight-Edge Diffraction Using Fresnel Approximation---------*/

extern complex float
Left_Straightedge_Diffraction_Float(float x, float edge, float F);

extern complex double
Left_Straightedge_Diffraction_Double(double x, double edge, double F);

extern complex long double
Left_Straightedge_Diffraction_Long_Double(long double x, long double edge,
                                          long double F);

extern complex float
Left_Straightedge_Diffraction_Float(float x, float edge, float F);

extern complex double
Left_Straightedge_Diffraction_Double(double x, double edge, double F);

extern complex long double
Left_Straightedge_Diffraction_Long_Double(long double x, long double edge,
                                          long double F);

extern double
Left_Straightedge_Diffraction_Char(char x, double a, double F);

extern double
Left_Straightedge_Diffraction_UChar(unsigned char x, double a, double F);

extern double
Left_Straightedge_Diffraction_Short(short x, double a, double F);

extern double
Left_Straightedge_Diffraction_UShort(unsigned short x, double a, double F);

extern double
Left_Straightedge_Diffraction_Int(int x, double a, double F);

extern double
Left_Straightedge_Diffraction_UInt(unsigned int x, double a, double F);

extern double
Left_Straightedge_Diffraction_Long(long x, double a, double F);

extern double
Left_Straightedge_Diffraction_ULong(unsigned long x, double a, double F);

extern double
Left_Straightedge_Diffraction_Long_Long(long long x, double a, double F);

extern double
Left_Straightedge_Diffraction_ULong_Long(unsigned long long x,
                                         double a, double F);

/*------------Square Wave Diffraction Using Fresnel Approximation-------------*/
extern complex double
Square_Wave_Diffraction_Double(double x, double W, double F, long N);
#endif