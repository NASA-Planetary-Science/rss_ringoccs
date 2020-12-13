/*  Include guard to avoid importing this file twice.                         */
#ifndef __RSS_RINGOCCS_RECONSTRUCTION_H__
#define __RSS_RINGOCCS_RECONSTRUCTION_H__

/*  Various functions, complex variables, and more found here.                */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>

/*  Structure that contains all of the necessary data.                        */
typedef struct TAUOBj {
    rssringoccs_ComplexDouble *T_in;
    double *rho_km_vals;
    double *F_km_vals;
    double *phi_rad_vals;
    double *kd_vals;
    double *B_rad_vals;
    double *D_km_vals;
    double *w_km_vals;
    double *t_oet_spm_vals;
    double *t_ret_spm_vals;
    double *t_set_spm_vals;
    double *rho_corr_pole_km_vals;
    double *rho_corr_timing_km_vals;
    double *tau_threshold_vals;
    double *phi_rl_rad_vals;
    double *p_norm_vals;
    double *phase_rad_vals;
    double ecc;
    double peri;
    double perturb[5];
    unsigned long start;
    unsigned long n_used;
    unsigned long arr_size;
    rssringoccs_Bool use_norm;
    rssringoccs_Bool use_fwd;
    rssringoccs_Bool use_fft;
    rssringoccs_Bool error_occurred;
    const char *error_message;
    const char *wtype;
    const char *psitype;
    unsigned char order;
    unsigned char interp;
    unsigned char status;
    rssringoccs_ComplexDouble *T_out;
} TAUObj;

/******************************************************************************
 *  Functions defined in _diffraction_functions.c                             *
 ******************************************************************************/
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
                                            double *coeffs, double dx, double F,
                                            double kd, unsigned long n_pts,
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
                                           double *coeffs, double dx, double F,
                                           double kd, unsigned long n_pts,
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
                                     unsigned long toler, double dx,
                                     double F, unsigned long n_pts,
                                     unsigned long center);

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
                                        unsigned long toler, double dx,
                                        double F, unsigned long n_pts,
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
                                    unsigned long toler, double dx, double F,
                                    unsigned long n_pts, unsigned long center);

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
                                      unsigned long toler, double dx,
                                      double F, unsigned long n_pts,
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
                                               double dx, double F,
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
                                      unsigned long toler, double dx,
                                      double F, unsigned long n_pts,
                                      unsigned long center, double ecc,
                                      double peri);

/*  Functions that compute the Fresnel Transform on a DLPObj instance.        */
extern void DiffractionCorrectionFresnel(TAUObj *dlp);
extern void DiffractionCorrectionLegendre(TAUObj *dlp);
extern void DiffractionCorrectionNewton(TAUObj *dlp);
extern void DiffractionCorrectionPerturbedNewton(TAUObj *dlp);
extern void DiffractionCorrectionEllipse(TAUObj *dlp);
extern void DiffractionCorrectionSimpleFFT(TAUObj *dlp);

#endif
