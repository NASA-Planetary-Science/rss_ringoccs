/*  Include guard to avoid importing this file twice.                         */
#ifndef RSS_RINGOCCS_RECONSTRUCTION_H
#define RSS_RINGOCCS_RECONSTRUCTION_H

/*  Various functions, complex variables, and more found here.                */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_history.h>

/*  size_t typedef is given here.                                             */
#include <stdlib.h>

typedef double (*rssringoccs_window_func)(double, double);

typedef enum {
    rssringoccs_DR_Fresnel,
    rssringoccs_DR_Legendre,
    rssringoccs_DR_Newton,
    rssringoccs_DR_NewtonD,
    rssringoccs_DR_NewtonDOld,
    rssringoccs_DR_NewtonDPhi,
    rssringoccs_DR_NewtonPerturb,
    rssringoccs_DR_Quadratic,
    rssringoccs_DR_Cubic,
    rssringoccs_DR_Quartic,
    rssringoccs_DR_QuarticD,
    rssringoccs_DR_SimpleFFT,
    rssringoccs_DR_Elliptical,
    rssringoccs_DR_None
} rssringoccs_Psitype_Enum;

/*  Structure that contains all of the necessary data.                        */
typedef struct rssringoccs_TAUObj_Def {
    tmpl_ComplexDouble *T_in;
    tmpl_ComplexDouble *T_out;
    tmpl_ComplexDouble *T_fwd;
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
    double *p_norm_fwd_vals;
    double *power_vals;
    double *phase_rad_vals;
    double *phase_fwd_vals;
    double *phase_vals;
    double *tau_fwd_vals;
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
    unsigned int toler;
    size_t start;
    size_t n_used;
    size_t arr_size;
    rssringoccs_window_func window_func;
    rssringoccs_Psitype_Enum psinum;
    tmpl_Bool use_norm;
    tmpl_Bool use_fwd;
    tmpl_Bool bfac;
    tmpl_Bool verbose;
    tmpl_Bool error_occurred;
    char *error_message;
    char *wtype;
    char *psitype;
    unsigned int order;
    rssringoccs_HistoryObj *history;
} rssringoccs_TAUObj;

typedef void (*rssringoccs_FresT)(rssringoccs_TAUObj *, double *,
                                  size_t, size_t);

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
                             double width, size_t nw_pts,
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

extern void
rssringoccs_Write_TAU_History(rssringoccs_TAUObj *tau);

#endif
