/*  Include guard to avoid importing this file twice.                         */
#ifndef __RSS_RINGOCCS_RECONSTRUCTION_H__
#define __RSS_RINGOCCS_RECONSTRUCTION_H__

/*  Various functions, complex variables, and more found here.                */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Structure that contains all of the necessary data.                        */
typedef struct TAUOBj {
    rssringoccs_ComplexDouble *T_in;
    double *rho_km_vals;
    double *F_km_vals;
    double *phi_rad_vals;
    double *kd_vals;
    double *f_sky_hz_vals;
    double *rho_dot_kms_vals;
    double *raw_tau_threshold_vals;
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
    char *error_message;
    char *wtype;
    char *psitype;
    unsigned char order;
    unsigned char interp;
    rssringoccs_ComplexDouble *T_out;
} rssringoccs_TAUObj;

extern void GetRangeFromString(char *range, double *rng_list);
extern void GetNormeqFromString(char *wtype, double *norm_eq);
extern void check_tau_data(rssringoccs_TAUObj *tau);
extern void check_tau_data_range(rssringoccs_TAUObj *dlp, double two_dx);

extern void select_window_func(rss_ringoccs_window_func *fw,
                               rssringoccs_TAUObj *tau);

extern void rssringoccs_Set_Tau_Error_Message(const char *mess,
                                              rssringoccs_TAUObj *tau);


extern void reset_window(double *x_arr, double *w_func, double dx, double width,
                         long nw_pts, rss_ringoccs_window_func fw);

/*  Functions that compute the Fresnel Transform on a TAUObj instance.        */
extern void DiffractionCorrectionFresnel(rssringoccs_TAUObj *tau);
extern void DiffractionCorrectionLegendre(rssringoccs_TAUObj *tau);
extern void DiffractionCorrectionNewton(rssringoccs_TAUObj *tau);
extern void DiffractionCorrectionEllipse(rssringoccs_TAUObj *tau);
extern void DiffractionCorrectionSimpleFFT(rssringoccs_TAUObj *tau);
extern void DiffractionCorrectionPerturbedNewton(rssringoccs_TAUObj *tau);

#endif
