/*  Include guard to avoid importing this file twice.                         */
#ifndef __RSS_RINGOCCS_RECONSTRUCTION_H__
#define __RSS_RINGOCCS_RECONSTRUCTION_H__

/*  Various functions, complex variables, and more found here.                */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

typedef enum {
    rssringoccs_DR_Fresnel,
    rssringoccs_DR_Legendre,
    rssringoccs_DR_Newton,
    rssringoccs_DR_NewtonD,
    rssringoccs_DR_NewtonDOld,
    rssringoccs_DR_NewtonDPhi,
    rssringoccs_DR_NewtonPerturb,
    rssringoccs_DR_Elliptical,
    rssringoccs_DR_None
} rssringoccs_Psitype_Enum;

/*  Structure that contains all of the necessary data.                        */
typedef struct rssringoccs_TAUObj {
    rssringoccs_ComplexDouble *T_in;
    rssringoccs_ComplexDouble *T_out;
    double *rho_km_vals;
    double *F_km_vals;
    double *phi_rad_vals;
    double *k_vals;
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
    double *power_vals;
    double *phase_rad_vals;
    double *phase_vals;
    double *tau_vals;
    double *rx_km_vals;
    double *ry_km_vals;
    double *rz_km_vals;
    double dx_km;
    double normeq;
    double sigma;
    double ecc;
    double peri;
    double res;
    double perturb[5];
    double rng_list[2];
    double rng_req[2];
    double EPS;
    unsigned char toler;
    unsigned long start;
    unsigned long n_used;
    unsigned long arr_size;
    rssringoccs_window_func window_func;
    rssringoccs_Psitype_Enum psinum;
    rssringoccs_Bool use_norm;
    rssringoccs_Bool use_fwd;
    rssringoccs_Bool use_fft;
    rssringoccs_Bool bfac;
    rssringoccs_Bool verbose;
    rssringoccs_Bool error_occurred;
    char *error_message;
    char *wtype;
    char *psitype;
    unsigned char order;
    unsigned char interp;
} rssringoccs_TAUObj;

typedef rssringoccs_ComplexDouble (*rssringoccs_FresT)(rssringoccs_TAUObj *,
                                                       double *,
                                                       unsigned long,
                                                       unsigned long);

extern void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Set_WType(const char *wtype, rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Set_Psitype(const char *psitype, rssringoccs_TAUObj* tau);

extern void
rssringoccs_Tau_Set_Range_From_String(const char *range,
                                      rssringoccs_TAUObj* tau);

extern rssringoccs_TAUObj *
rssringoccs_Create_TAUObj(rssringoccs_DLPObj *dlp, double res);

extern void
rssringoccs_Copy_DLP_Data_To_Tau(rssringoccs_DLPObj *dlp,
                                 rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Data_Range(rssringoccs_TAUObj *dlp);

extern void
rssringoccs_Tau_Check_Data(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Keywords(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Check_Occ_Type(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Compute_Vars(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Tau_Get_Window_Width(rssringoccs_TAUObj* tau);

extern void
rssringoccs_Tau_Finish(rssringoccs_TAUObj* tau);

extern void
rssringoccs_Destroy_Tau_Members(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Destroy_Tau(rssringoccs_TAUObj **tau);

extern void
rssringoccs_Tau_Reset_Window(double *x_arr, double *w_func, double dx,
                             double width, long nw_pts,
                             rssringoccs_window_func fw);

/*  Functions that compute the Fresnel Transform on a TAUObj instance.        */
extern void
rssringoccs_Diffraction_Correction_Fresnel(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_Legendre(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_Newton(rssringoccs_TAUObj *tau);

extern void
rssringoccs_Diffraction_Correction_SimpleFFT(rssringoccs_TAUObj *tau);

#endif
