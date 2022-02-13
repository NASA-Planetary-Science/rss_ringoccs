/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN

/*  The following are NON-STANDARD header files for using the C-Python API,   *
 *  and the numpy API.                                                        */
#include <Python.h>
#include <structmember.h>

/*  The standard library header stdlib contains malloc, calloc, and realloc.  */
#include <stdlib.h>

/*  The puts function is found here.                                          */
#include <stdio.h>

/*  The following header files come from libtmpl. This is a NON-STANDARD      *
 *  library that contains a variety of math tools.                            */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/tmpl_special_functions.h>

/*  And a bunch of headers from this project.                                 */
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  The definition of the DiffractionCorrection class as a C struct.          */
typedef struct PyDiffrecObj_Def {
    PyObject_HEAD
    PyObject *B_rad_vals;
    PyObject *D_km_vals;
    PyObject *F_km_vals;
    PyObject *f_sky_hz_vals;
    PyObject *p_norm_fwd_vals;
    PyObject *p_norm_vals;
    PyObject *phase_fwd_vals;
    PyObject *phase_rad_vals;
    PyObject *phase_vals;
    PyObject *phi_rad_vals;
    PyObject *phi_rl_rad_vals;
    PyObject *power_vals;
    PyObject *raw_tau_threshold_vals;
    PyObject *rev_info;
    PyObject *input_vars;
    PyObject *input_kwds;
    PyObject *rho_corr_pole_km_vals;
    PyObject *rho_corr_timing_km_vals;
    PyObject *rho_dot_kms_vals;
    PyObject *rho_km_vals;
    PyObject *t_oet_spm_vals;
    PyObject *t_ret_spm_vals;
    PyObject *t_set_spm_vals;
    PyObject *tau_threshold_vals;
    PyObject *tau_vals;
    PyObject *tau_fwd_vals;
    PyObject *w_km_vals;
    PyObject *history;
    PyObject *rx_km_vals;
    PyObject *ry_km_vals;
    PyObject *rz_km_vals;
    tmpl_Bool bfac;
    tmpl_Bool use_fwd;
    tmpl_Bool use_norm;
    tmpl_Bool verbose;
    double ecc;
    double input_res;
    double peri;
    double res_factor;
    double sigma;
    const char *psitype;
    const char *wtype;
} PyDiffrecObj;

/*  The definition of the DiffractionCorrection class as a C struct.          */
typedef struct PyCSVObj_Def {
    PyObject_HEAD
    PyObject *B_rad_vals;
    PyObject *D_km_vals;
    PyObject *tau_vals;
    PyObject *f_sky_hz_vals;
    PyObject *p_norm_vals;
    PyObject *phase_rad_vals;
    PyObject *phi_rad_vals;
    PyObject *phi_rl_rad_vals;
    PyObject *raw_tau_threshold_vals;
    PyObject *rev_info;
    PyObject *input_vars;
    PyObject *input_kwds;
    PyObject *rho_corr_pole_km_vals;
    PyObject *rho_corr_timing_km_vals;
    PyObject *rho_dot_kms_vals;
    PyObject *rho_km_vals;
    PyObject *t_oet_spm_vals;
    PyObject *t_ret_spm_vals;
    PyObject *t_set_spm_vals;
    PyObject *history;
    PyObject *rx_km_vals;
    PyObject *ry_km_vals;
    PyObject *rz_km_vals;
    PyObject *tau_rho;
    PyObject *tau_phase;
    PyObject *tau_power;
} PyCSVObj;


extern void
set_var(PyObject **py_ptr, double **ptr, unsigned long int len);

extern void capsule_cleanup(PyObject *capsule);

extern double *
extract_data(rssringoccs_DLPObj *dlp, PyObject *py_dlp, const char *var_name);

extern rssringoccs_DLPObj *
rssringoccs_Py_DLP_to_C_DLP(PyObject *py_dlp);

extern void
rssringoccs_C_Tau_to_Py_Tau(PyDiffrecObj *py_tau, rssringoccs_TAUObj *tau);

extern void
rssringoccs_Get_Py_Perturb(rssringoccs_TAUObj *tau, PyObject *perturb);

extern void rssringoccs_Get_Py_Range(rssringoccs_TAUObj *tau, PyObject *rngreq);

extern void
rssringoccs_Get_Py_Vars_From_Tau_Self(rssringoccs_TAUObj *tau,
                                      PyDiffrecObj *self);

extern void Diffrec_dealloc(PyDiffrecObj *self);

extern int Diffrec_init(PyDiffrecObj *self, PyObject *args, PyObject *kwds);

extern PyTypeObject DiffrecType;

extern void
rssringoccs_C_CSV_to_Py_CSV(PyCSVObj *py_csv, rssringoccs_CSVData *csv);

extern void ExtractCSVData_dealloc(PyCSVObj *self);

extern int
ExtractCSVData_init(PyCSVObj *self, PyObject *args, PyObject *kwds);

extern PyTypeObject ExtractCSVDataType;

