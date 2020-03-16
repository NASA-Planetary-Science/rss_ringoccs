/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_KERNEL_H
#define RSS_RINGOCCS_FRESNEL_KERNEL_H

/*  Various mathematical functions defined here.                              */
#include <math.h>

/*---------------------------Newton-Raphson Function--------------------------*/

extern float
Newton_Raphson_Float(float x, float (*f)(float), float (*f_prime)(float),
                     float EPS, long toler);

extern double
Newton_Raphson_Double(double x, double (*f)(double),
                      double (*f_prime)(double), double EPS, long toler);

extern long double
Newton_Raphson_Long_Double(long double x, long double (*f)(long double),
                           long double (*f_prime)(long double), long double EPS,
                           long toler);

/*----------------------------The Fresnel Kernel------------------------------*/

extern float Fresnel_Psi_Float(float kD, float r, float r0, float phi,
                               float phi0, float B, float D);

extern double Fresnel_Psi_Double(double kD, double r, double r0, double phi,
                                 double phi0, double B, double D);

extern double Fresnel_Psi_Long_Double(long double kD, long double r,
                                      long double r0, long double phi,
                                      long double phi0, long double B,
                                      long double D);

/*----------------The First Derivative of the Fresnel Kernel------------------*/

extern float Fresnel_dPsi_dPhi_Float(float kD, float r, float r0, float phi,
                                     float phi0, float B, float D);

extern double Fresnel_dPsi_dPhi_Double(double kD, double r, double r0,
                                       double phi, double phi0, double B,
                                       double D);

extern double Fresnel_dPsi_dPhi_Long_Double(long double kD, long double r,
                                            long double r0, long double phi,
                                            long double phi0, long double B,
                                            long double D);

/*------The Derivative of the Fresnel Kernel with Elliptic Perturbations------*/

extern float Fresnel_dPsi_dPhi_Ellipse_Float(float kD,  float r,    float r0,
                                             float phi, float phi0, float B,
                                             float D,   float ecc,  float peri);

extern double
Fresnel_dPsi_dPhi_Ellipse_Double(double kD,   double r, double r0, double phi,
                                 double phi0, double B, double D,  double ecc,
                                 double peri);

extern long double
Fresnel_dPsi_dPhi_Ellipse_Long_Double(long double kD,   long double r,
                                      long double r0,   long double phi,
                                      long double phi0, long double B,
                                      long double D,    long double ecc,
                                      long double peri);

/*---------------The Second Derivative of the Fresnel Kernel------------------*/

extern float Fresnel_d2Psi_dPhi2_Float(float kD,   float r, float r0, float phi,
                                       float phi0, float B, float D);

extern double
Fresnel_d2Psi_dPhi2_Double(double kD,  double r,    double r0,
                           double phi, double phi0, double B, double D);

extern long double
Fresnel_d2Psi_dPhi2_Long_Double(long double kD, long double r,
                                long double r0, long double phi,
                                long double phi0, long double B, long double D);

/*--------------Functions to Perform Newton-Raphon on Psi---------------------*/

extern double Newton_Raphson_Fresnel_Psi(double kD,  double r,    double r0,
                                         double phi, double phi0, double B,
                                         double D,   double EPS, long toler);

extern double
Newton_Raphson_Fresnel_Ellipse(double kD,   double r, double r0, double phi,
                               double phi0, double B, double D,  double ecc,
                               double peri, double EPS, unsigned char toler);

#endif
