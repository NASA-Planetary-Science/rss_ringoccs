/*  Include guard to avoid importing this file twice.                         */
#ifndef _RSS_RINGOCCS_FRESNEL_KERNEL_H_
#define _RSS_RINGOCCS_FRESNEL_KERNEL_H_

/*----------------------------The Fresnel Kernel------------------------------*/
extern float
rssringoccs_Float_Fresnel_Psi(float k, float r, float r0, float phi,
                              float phi0, float B, float D);

extern double
rssringoccs_Double_Fresnel_Psi(double k, double r, double r0, double phi,
                               double phi0, double B, double D);

extern long double
rssringoccs_LDouble_Fresnel_Psi(long double k, long double r,
                                long double r0, long double phi,
                                long double phi0, long double B, long double D);

extern double rssringoccs_Double_Fresnel_Psi_Old(double kD, double r, double r0,
                                                 double phi, double phi0,
                                                 double B, double D);

extern double Newton_Raphson_Fresnel_Psi_D(double k, double r, double r0,
                                           double phi, double phi0, double B,
                                           double EPS, long toler, double rx,
                                           double ry, double rz);

extern double Newton_Raphson_Fresnel_Psi_dD_dphi(double k, double r, double r0,
                                                 double phi, double phi0,
                                                 double B, double EPS,
                                                 long toler, double rx,
                                                 double ry, double rz);

extern double
Newton_Raphson_Fresnel_Psi_D_Old(double kD, double r, double r0, double phi,
                                 double phi0, double B, double EPS, long toler,
                                 double rx, double ry, double rz);

extern double
rssringoccs_Double_Fresnel_dPsi_dPhi_D(double k, double r, double r0,
                                       double phi, double phi0, double B,
                                       double rx, double ry, double rz);

/*----------------The First Derivative of the Fresnel Kernel------------------*/
extern float
rssringoccs_Float_Fresnel_dPsi_dPhi(float k, float r, float r0, float phi,
                                    float phi0, float B, float D);

extern double
rssringoccs_Double_Fresnel_dPsi_dPhi(double k, double r, double r0, double phi,
                                     double phi0, double B, double D);

extern long double
rssringoccs_LDouble_Fresnel_dPsi_dPhi(long double k, long double r,
                                      long double r0, long double phi,
                                      long double phi0, long double B,
                                      long double D);

/*------The Derivative of the Fresnel Kernel with Elliptic Perturbations------*/
extern float
rssringoccs_Float_Fresnel_dPsi_dPhi_Ellipse(float k,  float r,    float r0,
                                            float phi, float phi0, float B,
                                            float D,   float ecc,  float peri);

extern double
rssringoccs_Double_Fresnel_dPsi_dPhi_Ellipse(double k, double r, double r0,
                                             double phi, double phi0, double B,
                                             double D, double ecc, double peri);

extern long double
rssringoccs_LDouble_Fresnel_dPsi_dPhi_Ellipse(long double k,   long double r,
                                              long double r0,   long double phi,
                                              long double phi0, long double B,
                                              long double D,    long double ecc,
                                              long double peri);

/*---------------The Second Derivative of the Fresnel Kernel------------------*/
extern float
rssringoccs_Float_Fresnel_d2Psi_dPhi2(float k,   float r, float r0, float phi,
                                      float phi0, float B, float D);

extern double
rssringoccs_Double_Fresnel_d2Psi_dPhi2(double k,  double r,    double r0,
                                       double phi, double phi0, double B,
                                       double D);

extern long double
rssringoccs_LDouble_Fresnel_d2Psi_dPhi2(long double k, long double r,
                                        long double r0, long double phi,
                                        long double phi0, long double B,
                                        long double D);

/*--------------Functions to Perform Newton-Raphon on Psi---------------------*/
extern double
Newton_Raphson_Fresnel_Psi(double k,   double r, double r0, double phi,
                           double phi0, double B, double D,  double EPS,
                           unsigned char toler);

extern double
rssringoccs_Double_Newton_Raphson_Fresnel_Ellipse(double k, double r,
                                                  double r0, double phi,
                                                  double phi0, double B,
                                                  double ecc, double peri,
                                                  double EPS,
                                                  unsigned char toler,
                                                  double rx, double ry,
                                                  double rz);

/*  Fresnel scale.  */
extern float
rssringoccs_Float_Fresnel_Scale(float lambda, float d, float phi, float b);

extern double
rssringoccs_Double_Fresnel_Scale(double lambda, double d, double phi, double b);

extern long double
rssringoccs_LDouble_Fresnel_Scale(long double lambda, long double d,
                                  long double phi, long double b);

#endif
