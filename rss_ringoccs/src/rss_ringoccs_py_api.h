#ifndef RSS_RINGOCCS_H
#define RSS_RINGOCCS_H

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>

#include <structmember.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>

typedef struct rssringoccs_Generic_Function_Obj {
    long (*long_func)(long);
    float (*float_func)(float);
    double (*double_func)(double);
    long double (*ldouble_func)(long double);
    rssringoccs_ComplexDouble (*cdouble_from_real_func)(double);
    rssringoccs_ComplexDouble
        (*cdouble_from_complex_func)(rssringoccs_ComplexDouble);
    const char *func_name;
} rssringoccs_Generic_Function_Obj;

typedef struct {
    PyObject_HEAD
    PyObject          *B_rad_vals;
    PyObject          *D_km_vals;
    PyObject          *F_km_vals;
    PyObject          *f_sky_hz_vals;
    PyObject          *p_norm_fwd_vals;
    PyObject          *p_norm_vals;
    PyObject          *phase_fwd_vals;
    PyObject          *phase_rad_vals;
    PyObject          *phase_vals;
    PyObject          *phi_rad_vals;
    PyObject          *phi_rl_rad_vals;
    PyObject          *power_vals;
    PyObject          *T_hat_vals;
    PyObject          *T_hat_fwd_vals;
    PyObject          *T_vals;
    PyObject          *raw_tau_threshold_vals;
    PyObject          *rev_info;
    PyObject          *rho_corr_pole_km_vals;
    PyObject          *rho_corr_timing_km_vals;
    PyObject          *rho_dot_kms_vals;
    PyObject          *rho_km_vals;
    PyObject          *t_oet_spm_vals;
    PyObject          *t_ret_spm_vals;
    PyObject          *t_set_spm_vals;
    PyObject          *tau_threshold_vals;
    PyObject          *tau_vals;
    PyObject          *w_km_vals;
    PyObject          *dathist;
    PyObject          *history;
    PyObject          *rx_km_vals;
    PyObject          *ry_km_vals;
    PyObject          *rz_km_vals;
    rssringoccs_Bool   bfac;
    rssringoccs_Bool   use_fwd;
    rssringoccs_Bool   use_norm;
    rssringoccs_Bool   verbose;
    rssringoccs_Bool   write_file;
    double             ecc;
    double             input_res;
    double             peri;
    double             res_factor;
    double             sigma;
    const char        *psitype;
    const char        *wtype;
    unsigned char      interp;
} PyDiffrecObj;

#endif
