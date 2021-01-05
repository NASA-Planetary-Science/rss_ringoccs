/*  Include guard to avoid importing this file twice.                         */
#ifndef __RSS_RINGOCCS_FRESNEL_TRANFORM_H__
#define __RSS_RINGOCCS_FRESNEL_TRANFORM_H__

/*  Various functions, complex variables, and more found here.                */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

extern rssringoccs_ComplexDouble
Fresnel_Transform_Double(double *x_arr, rssringoccs_ComplexDouble *T_in,
                         double *w_func, double F, double dx,
                         unsigned long n_pts, unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Norm_Double(double *x_arr, rssringoccs_ComplexDouble *T_in,
                              double *w_func, double F,
                              unsigned long n_pts, unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Legendre_Even_Double(double *x_arr,
                                       rssringoccs_ComplexDouble *T_in,
                                       double *w_func, double D, double *coeffs,
                                       double dx, double F, double kd,
                                       unsigned long n_pts, unsigned char order,
                                       unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Legendre_Norm_Even_Double(double *x_arr,
                                            rssringoccs_ComplexDouble *T_in,
                                            double *w_func, double D,
                                            double *coeffs, double kd,
                                            unsigned long n_pts,
                                            unsigned char order,
                                            unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Legendre_Odd_Double(double *x_arr,
                                      rssringoccs_ComplexDouble *T_in,
                                      double *w_func, double D, double *coeffs,
                                      double dx, double F, double kd,
                                      unsigned long n_pts, unsigned char order,
                                      unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Legendre_Norm_Odd_Double(double *x_arr,
                                           rssringoccs_ComplexDouble *T_in,
                                           double *w_func, double D,
                                           double *coeffs, double kd,
                                           unsigned long n_pts,
                                           unsigned char order,
                                           unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Newton_Double(double *x_arr, double *phi_arr,
                                rssringoccs_ComplexDouble *T_in, double *w_func,
                                double kD, double r, double B, double D,
                                double EPS, unsigned long toler, double dx,
                                double F, unsigned long n_pts,
                                unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Newton_Norm_Double(double *x_arr, double *phi_arr,
                                     rssringoccs_ComplexDouble *T_in,
                                     double *w_func, double kD, double r,
                                     double B, double D, double EPS,
                                     unsigned long toler, unsigned long n_pts,
                                     unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Newton_D_Double(double *x_arr, double *phi_arr,
                                  rssringoccs_ComplexDouble *T_in, double *w_func,
                                  double kD, double r, double B, double EPS,
                                  unsigned long toler, double dx,
                                  double F, unsigned long n_pts,
                                  unsigned long center, double rx, double ry,
                                  double rz);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Newton_D_Norm_Double(double *x_arr, double *phi_arr,
                                       rssringoccs_ComplexDouble *T_in,
                                       double *w_func, double k, double r,
                                       double B, double EPS, long toler,
                                       long n_pts, long center,
                                       double rx, double ry, double rz);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Newton_D_Old_Double(double *x_arr, double *phi_arr,
                                  rssringoccs_ComplexDouble *T_in, double *w_func,
                                  double kD, double r, double B, double EPS,
                                  unsigned long toler, double dx,
                                  double F, unsigned long n_pts,
                                  unsigned long center, double rx, double ry,
                                  double rz);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Newton_D_Old_Norm_Double(double *x_arr, double *phi_arr,
                                       rssringoccs_ComplexDouble *T_in,
                                       double *w_func, double k, double r,
                                       double B, double EPS, long toler,
                                       long n_pts, long center,
                                       double rx, double ry, double rz);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Quadratic_Double(double *x_arr, double *phi_arr,
                                   rssringoccs_ComplexDouble *T_in,
                                   double *w_func, double kD, double r,
                                   double B, double D, double EPS,
                                   unsigned long toler, double dx, double F,
                                   unsigned long n_pts, unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Quadratic_Norm_Double(double *x_arr, double *phi_arr,
                                        rssringoccs_ComplexDouble *T_in,
                                        double *w_func, double kD, double r,
                                        double B, double D, double EPS,
                                        unsigned long toler,
                                        unsigned long n_pts,
                                        unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Cubic_Double(double *x_arr, double *phi_arr,
                               rssringoccs_ComplexDouble *T_in, double *w_func,
                               double kD, double r, double B, double D,
                               double EPS, unsigned long toler, double dx,
                               double F, unsigned long n_pts,
                               unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Cubic_Norm_Double(double *x_arr, double *phi_arr,
                                    rssringoccs_ComplexDouble *T_in,
                                    double *w_func, double kD, double r,
                                    double B, double D, double EPS,
                                    unsigned long toler,
                                    unsigned long n_pts,
                                    unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Quartic_Double(double *x_arr, double *phi_arr,
                                 rssringoccs_ComplexDouble *T_in,
                                 double *w_func, double kD, double r, double B,
                                 double D, double EPS, unsigned long toler,
                                 double dx, double F, unsigned long n_pts,
                                 unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Quartic_Norm_Double(double *x_arr, double *phi_arr,
                                      rssringoccs_ComplexDouble *T_in,
                                      double *w_func, double kD, double r,
                                      double B, double D, double EPS,
                                      unsigned long toler, unsigned long n_pts,
                                      unsigned long center);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Perturbed_Newton_Double(double *x_arr, double *phi_arr,
                                          rssringoccs_ComplexDouble *T_in,
                                          double *w_func, double kD, double r,
                                          double B, double D, double EPS,
                                          unsigned long toler, double dx,
                                          double F, unsigned long n_pts,
                                          unsigned long center,
                                          double perturb[5]);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Perturbed_Newton_Norm_Double(double *x_arr, double *phi_arr,
                                               rssringoccs_ComplexDouble *T_in,
                                               double *w_func, double kD,
                                               double r, double B, double D,
                                               double EPS, unsigned long toler,
                                               unsigned long n_pts,
                                               unsigned long center,
                                               double perturb[5]);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Ellipse_Double(double *x_arr, double *phi_arr,
                                 rssringoccs_ComplexDouble *T_in,
                                 double *w_func, double kD, double r, double B,
                                 double D, double EPS, unsigned long toler,
                                 double dx, double F, unsigned long n_pts,
                                 unsigned long center, double ecc, double peri);

extern rssringoccs_ComplexDouble
Fresnel_Transform_Ellipse_Norm_Double(double *x_arr, double *phi_arr,
                                      rssringoccs_ComplexDouble *T_in,
                                      double *w_func, double kD, double r,
                                      double B, double D, double EPS,
                                      unsigned long toler, unsigned long n_pts,
                                      unsigned long center, double ecc,
                                      double peri);

#endif
