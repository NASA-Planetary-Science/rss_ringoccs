/*  Include guard to avoid importing this file twice.                         */
#ifndef __RSS_RINGOCCS_FRESNEL_TRANFORM_H__
#define __RSS_RINGOCCS_FRESNEL_TRANFORM_H__

/*  Various functions, complex variables, and more found here.                */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

extern void
Fresnel_Transform_Double(rssringoccs_TAUObj *tau, double *x_arr, double *w_func,
                         unsigned long n_pts, unsigned long center);

extern void
Fresnel_Transform_Norm_Double(rssringoccs_TAUObj *tau, double *x_arr,
                              double *w_func, unsigned long n_pts,
                              unsigned long center);

extern void
Fresnel_Transform_Legendre_Even_Double(rssringoccs_TAUObj *tau, double *x_arr,
                                       double *w_func, double *coeffs,
                                       unsigned long n_pts,
                                       unsigned long center);

extern void
Fresnel_Transform_Legendre_Norm_Even_Double(rssringoccs_TAUObj *tau,
                                            double *x_arr, double *w_func,
                                            double *coeffs, unsigned long n_pts,
                                            unsigned long center);

extern void
Fresnel_Transform_Legendre_Odd_Double(rssringoccs_TAUObj *tau, double *x_arr,
                                      double *w_func, double *coeffs,
                                      unsigned long n_pts,
                                      unsigned long center);

extern void
Fresnel_Transform_Legendre_Norm_Odd_Double(rssringoccs_TAUObj *tau,
                                           double *x_arr, double *w_func,
                                           double *coeffs, unsigned long n_pts,
                                           unsigned long center);

extern void
Fresnel_Transform_Newton_Double(rssringoccs_TAUObj *tau,
                                double *w_func,
                                unsigned long n_pts,
                                unsigned long center);

extern void
Fresnel_Transform_Newton_Norm_Double(rssringoccs_TAUObj *tau,
                                     double *w_func,
                                     unsigned long n_pts,
                                     unsigned long center);

extern void
Fresnel_Transform_Newton_D_Double(rssringoccs_TAUObj *tau,
                                  double *w_func,
                                  unsigned long n_pts,
                                  unsigned long center);

extern void
Fresnel_Transform_Newton_D_Norm_Double(rssringoccs_TAUObj *tau,
                                       double *w_func,
                                       unsigned long n_pts,
                                       unsigned long center);

extern void
Fresnel_Transform_Newton_D_Old_Double(rssringoccs_TAUObj *tau,
                                      double *w_func,
                                      unsigned long n_pts,
                                      unsigned long center);

extern void
Fresnel_Transform_Newton_D_Old_Norm_Double(rssringoccs_TAUObj *tau,
                                           double *w_func,
                                           unsigned long n_pts,
                                           unsigned long center);

extern void
Fresnel_Transform_Newton_dD_dphi_Double(rssringoccs_TAUObj *tau,
                                        double *w_func,
                                        unsigned long n_pts,
                                        unsigned long center);

extern void
Fresnel_Transform_Newton_dD_dphi_Norm_Double(rssringoccs_TAUObj *tau,
                                             double *w_func,
                                             unsigned long n_pts,
                                             unsigned long center);

extern void
Fresnel_Transform_Perturbed_Newton_Double(rssringoccs_TAUObj *tau,
                                          double *w_func,
                                          unsigned long n_pts,
                                          unsigned long center);

extern void
Fresnel_Transform_Perturbed_Newton_Norm_Double(rssringoccs_TAUObj *tau,
                                               double *w_func,
                                               unsigned long n_pts,
                                               unsigned long center);

extern void
Fresnel_Transform_Quadratic_Double(rssringoccs_TAUObj *tau,
                                   double *w_func,
                                   unsigned long n_pts,
                                   unsigned long center);

extern void
Fresnel_Transform_Quadratic_Norm_Double(rssringoccs_TAUObj *tau,
                                        double *w_func,
                                        unsigned long n_pts,
                                        unsigned long center);

extern void
Fresnel_Transform_Cubic_Double(rssringoccs_TAUObj *tau,
                               double *w_func,
                               unsigned long n_pts,
                               unsigned long center);

extern void
Fresnel_Transform_Cubic_Norm_Double(rssringoccs_TAUObj *tau,
                                    double *w_func,
                                    unsigned long n_pts,
                                    unsigned long center);

extern void
Fresnel_Transform_Quartic_Double(rssringoccs_TAUObj *tau,
                                 double *w_func,
                                 unsigned long n_pts,
                                 unsigned long center);

extern void
Fresnel_Transform_Quartic_Norm_Double(rssringoccs_TAUObj *tau,
                                      double *w_func,
                                      unsigned long n_pts,
                                      unsigned long center);

extern void
Fresnel_Transform_Quartic_D_Double(rssringoccs_TAUObj *tau,
                                   double *w_func,
                                   unsigned long n_pts,
                                   unsigned long center);

extern void
Fresnel_Transform_Quartic_D_Norm_Double(rssringoccs_TAUObj *tau,
                                        double *w_func,
                                        unsigned long n_pts,
                                        unsigned long center);

extern void
Fresnel_Transform_Ellipse_Double(rssringoccs_TAUObj *tau,
                                 double *w_func,
                                 unsigned long n_pts,
                                 unsigned long center);

extern void
Fresnel_Transform_Ellipse_Norm_Double(rssringoccs_TAUObj *tau,
                                      double *w_func,
                                      unsigned long n_pts,
                                      unsigned long center);

#endif
