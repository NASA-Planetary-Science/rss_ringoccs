/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_KERNEL_H
#define RSS_RINGOCCS_FRESNEL_KERNEL_H

extern double
rssringoccs_Fresnel_Psi(double k, double r, double r0, double phi,
                        double phi0, double B, double D);

extern double
rssringoccs_Fresnel_Psi_Old(double kD, double r, double r0, double phi,
                            double phi0, double B, double D);

extern double
rssringoccs_Newton_Raphson_Fresnel_Psi_D(double k, double r, double r0,
                                         double phi, double phi0, double B,
                                         double EPS, unsigned int toler,
                                         double rx, double ry, double rz);

extern double
rssringoccs_Newton_Raphson_Fresnel_Psi_dD_dphi(double k, double r, double r0,
                                               double phi, double phi0,
                                               double B, double EPS,
                                               unsigned int toler, double rx,
                                               double ry, double rz);

extern double
rssringoccs_Newton_Raphson_Fresnel_Psi_D_Old(double kD, double r, double r0,
                                             double phi, double phi0, double B,
                                             double EPS, unsigned int toler,
                                             double rx, double ry, double rz);

extern double
rssringoccs_Fresnel_dPsi_dPhi_D(double k, double r, double r0, double phi,
                                double phi0, double B, double rx,
                                double ry, double rz);

extern double
rssringoccs_Fresnel_dPsi_dPhi(double k, double r, double r0, double phi,
                              double phi0, double B, double D);

extern double
rssringoccs_Fresnel_dPsi_dPhi_Ellipse(double k, double r, double r0, double phi,
                                      double phi0, double B, double D,
                                      double ecc, double peri);

extern double
rssringoccs_Fresnel_d2Psi_dPhi2(double k, double r, double r0, double phi,
                                double phi0, double B, double D);

extern double
rssringoccs_Newton_Raphson_Fresnel_Psi(double k, double r, double r0,
                                       double phi, double phi0, double B,
                                       double D,  double EPS,
                                       unsigned int toler);

extern double
rssringoccs_Newton_Raphson_Fresnel_Ellipse(double k, double r, double r0,
                                           double phi, double phi0, double B,
                                           double ecc, double peri, double EPS,
                                           unsigned int toler, double rx,
                                           double ry, double rz);

extern double
rssringoccs_Fresnel_Scale(double lambda, double d, double phi, double b);

#endif
