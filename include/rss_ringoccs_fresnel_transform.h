/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_FRESNEL_TRANFORM_H
#define RSS_RINGOCCS_FRESNEL_TRANFORM_H

/*  Various functions, complex variables, and more found here.                */
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  size_t typedef given here.                                                */
#include <stdlib.h>

extern void
rssringoccs_Fresnel_Transform(rssringoccs_TAUObj *tau, double *x_arr,
                              double *w_func, size_t n_pts, size_t center);

extern void
rssringoccs_Fresnel_Transform_Norm(rssringoccs_TAUObj *tau, double *x_arr,
                                   double *w_func, size_t n_pts, size_t center);

extern void
rssringoccs_Fresnel_Transform_Legendre_Even(rssringoccs_TAUObj *tau,
                                            double *x_arr, double *w_func,
                                            double *coeffs, size_t n_pts,
                                            size_t center);

extern void
rssringoccs_Fresnel_Transform_Legendre_Norm_Even(rssringoccs_TAUObj *tau,
                                                 double *x_arr, double *w_func,
                                                 double *coeffs, size_t n_pts,
                                                 size_t center);

extern void
rssringoccs_Fresnel_Transform_Legendre_Odd(rssringoccs_TAUObj *tau,
                                           double *x_arr, double *w_func,
                                           double *coeffs, size_t n_pts,
                                           size_t center);

extern void
rssringoccs_Fresnel_Transform_Legendre_Norm_Odd(rssringoccs_TAUObj *tau,
                                                double *x_arr, double *w_func,
                                                double *coeffs,
                                                size_t n_pts, size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton(rssringoccs_TAUObj *tau, double *w_func,
                                     size_t n_pts, size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton_Norm(rssringoccs_TAUObj *tau,
                                          double *w_func, size_t n_pts,
                                          size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton_D(rssringoccs_TAUObj *tau, double *w_func,
                                       size_t n_pts, size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton_D_Norm(rssringoccs_TAUObj *tau,
                                            double *w_func,
                                            size_t n_pts,
                                            size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton_D_Old(rssringoccs_TAUObj *tau,
                                           double *w_func,
                                           size_t n_pts,
                                           size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton_D_Old_Norm(rssringoccs_TAUObj *tau,
                                                double *w_func,
                                                size_t n_pts,
                                                size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton_dD_dphi(rssringoccs_TAUObj *tau,
                                             double *w_func,
                                             size_t n_pts,
                                             size_t center);

extern void
rssringoccs_Fresnel_Transform_Newton_dD_dphi_Norm(rssringoccs_TAUObj *tau,
                                                  double *w_func,
                                                  size_t n_pts,
                                                  size_t center);

extern void
rssringoccs_Fresnel_Transform_Perturbed_Newton(rssringoccs_TAUObj *tau,
                                               double *w_func,
                                               size_t n_pts,
                                               size_t center);

extern void
rssringoccs_Fresnel_Transform_Perturbed_Newton_Norm(rssringoccs_TAUObj *tau,
                                                    double *w_func,
                                                    size_t n_pts,
                                                    size_t center);

extern void
rssringoccs_Fresnel_Transform_Quadratic(rssringoccs_TAUObj *tau,
                                        double *w_func,
                                        size_t n_pts,
                                        size_t center);

extern void
rssringoccs_Fresnel_Transform_Quadratic_Norm(rssringoccs_TAUObj *tau,
                                             double *w_func,
                                             size_t n_pts,
                                             size_t center);

extern void
rssringoccs_Fresnel_Transform_Cubic(rssringoccs_TAUObj *tau,
                                    double *w_func,
                                    size_t n_pts,
                                    size_t center);

extern void
rssringoccs_Fresnel_Transform_Cubic_Norm(rssringoccs_TAUObj *tau,
                                         double *w_func,
                                         size_t n_pts,
                                         size_t center);

extern void
rssringoccs_Fresnel_Transform_Quartic(rssringoccs_TAUObj *tau,
                                      double *w_func,
                                      size_t n_pts,
                                      size_t center);

extern void
rssringoccs_Fresnel_Transform_Quartic_Norm(rssringoccs_TAUObj *tau,
                                           double *w_func,
                                           size_t n_pts,
                                           size_t center);

extern void
rssringoccs_Fresnel_Transform_Quartic_D(rssringoccs_TAUObj *tau,
                                        double *w_func,
                                        size_t n_pts,
                                        size_t center);

extern void
rssringoccs_Fresnel_Transform_Quartic_D_Norm(rssringoccs_TAUObj *tau,
                                             double *w_func,
                                             size_t n_pts,
                                             size_t center);

extern void
rssringoccs_Fresnel_Transform_Ellipse(rssringoccs_TAUObj *tau,
                                      double *w_func,
                                      size_t n_pts,
                                      size_t center);

extern void
rssringoccs_Fresnel_Transform_Ellipse_Norm(rssringoccs_TAUObj *tau,
                                           double *w_func,
                                           size_t n_pts,
                                           size_t center);

#endif
