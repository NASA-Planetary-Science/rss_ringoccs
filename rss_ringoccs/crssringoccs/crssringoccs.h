/******************************************************************************
 *                                  LICENSE                                   *
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

/*  The Python-C API is given here. The Python documentation recommends       *
 *  including Python.h before anything (even standard library headers).       */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

/*  Standard library file where size_t is declared.                           */
#include <stddef.h>

/*  tmpl_Bool typedef is given here. It provides Booleans for C89 compilers.  */
#include <libtmpl/include/tmpl_bool.h>

/*  And a bunch of headers from this project.                                 */
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  The definition of the DiffractionCorrection class as a C struct.          */
typedef struct PyDiffrecObj_Def {
    PyObject_HEAD
    PyObject *T_in;                   /*  Input complex transmittance.        */
    PyObject *T_out;                  /*  Reconstructed complex transmittance.*/
    PyObject *T_fwd;                  /*  Forward model complex transmittance.*/
    PyObject *k_vals;                 /*  Wavenumber, 2 pi / wavelength.      */
    PyObject *B_deg_vals;             /*  Ring opening angle.                 */
    PyObject *D_km_vals;              /*  Spacecraft-to-Ring distance.        */
    PyObject *F_km_vals;              /*  Fresnel scale.                      */
    PyObject *phi_deg_vals;           /*  Ring azimuth angle.                 */
    PyObject *phi_rl_deg_vals;        /*  Ring longitude angle.               */
    PyObject *input_vars;             /*  Input parameters for the class.     */
    PyObject *input_kwds;             /*  Input keywords for the class.       */
    PyObject *rho_corr_pole_km_vals;  /*  Pole corrected ring radius.         */
    PyObject *rho_corr_timing_km_vals;/*  Timing corrected ring radius.       */
    PyObject *rho_dot_kms_vals;       /*  Ring radial velocity.               */
    PyObject *rho_km_vals;            /*  Ring radius.                        */
    PyObject *t_oet_spm_vals;         /*  Seconds past midnight, observer.    */
    PyObject *t_ret_spm_vals;         /*  Seconds past midnight, ring.        */
    PyObject *t_set_spm_vals;         /*  Seconds past midngith, spacecraft.  */
    PyObject *tau_threshold_vals;     /*  Reconstructed tau threshold.        */
    PyObject *w_km_vals;              /*  Window width.                       */
    PyObject *rx_km_vals;             /*  x component of spacecraft.          */
    PyObject *ry_km_vals;             /*  y component of spacecraft.          */
    PyObject *rz_km_vals;             /*  z component of spacecraft.          */
    tmpl_Bool bfac;                   /*  Boolean for b factor in resolution. */
    tmpl_Bool use_fwd;                /*  Boolean for forward modeling.       */
    tmpl_Bool use_norm;               /*  Boolean for window normalization.   */
    tmpl_Bool verbose;                /*  Boolean for printing messages.      */
    double ecc;                       /*  Eccentricity, elliptical rings only.*/
    double input_res;                 /*  Input resolution, in kilometers.    */
    double peri;                      /*  Periapse, elliptical rings only.    */
    double res_factor;                /*  Resolution scale factor, unitless.  */
    double sigma;                     /*  Allen deviation of spacecraft.      */
    const char *outfiles;             /*  TAB files for this Tau object.      */
    const char *wtype;
    const char *psitype;
} PyDiffrecObj;

/*  The CSV struct containing all of the data for diffraction reconstruction. */
typedef struct PyCSVObj_Def {
    PyObject_HEAD
    PyObject *B_deg_vals;
    PyObject *D_km_vals;
    PyObject *f_sky_hz_vals;
    PyObject *p_norm_vals;
    PyObject *raw_tau_vals;
    PyObject *phase_deg_vals;
    PyObject *phi_deg_vals;
    PyObject *phi_rl_deg_vals;
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
    PyObject *tau_phase;
    PyObject *tau_power;
    PyObject *tau_vals;
} PyCSVObj;

extern void crssringoccs_Set_Var(PyObject **py_ptr,
                                 double *ptr,
                                 size_t len);

extern void crssringoccs_Set_CVar(PyObject **py_ptr,
                                  tmpl_ComplexDouble *ptr,
                                  size_t len);

extern void crssringoccs_Capsule_Cleanup(PyObject *capsule);

extern double *
crssringoccs_Extract_Data(rssringoccs_DLPObj *dlp,
                          PyObject *py_dlp,
                          const char *var_name);

extern rssringoccs_DLPObj *crssringoccs_Py_DLP_To_C_DLP(PyObject *py_dlp);

extern void
crssringoccs_C_Tau_To_Py_Tau(PyDiffrecObj *py_tau, rssringoccs_TAUObj *tau);

extern void
crssringoccs_Get_Py_Perturb(rssringoccs_TAUObj *tau, PyObject *perturb);

extern void
crssringoccs_Get_Py_Range(rssringoccs_TAUObj *tau, PyObject *rngreq);

extern void
crssringoccs_Get_Py_Vars_From_Tau_Self(rssringoccs_TAUObj *tau,
                                       const PyDiffrecObj *self);

extern void Diffrec_dealloc(PyDiffrecObj *self);

extern int Diffrec_init(PyDiffrecObj *self, PyObject *args, PyObject *kwds);

extern PyTypeObject DiffrecType;







extern void
crssringoccs_C_CSV_to_Py_CSV(PyCSVObj *py_csv, rssringoccs_CSVData *csv);

extern void ExtractCSVData_dealloc(PyCSVObj *self);

extern int
ExtractCSVData_init(PyCSVObj *self, PyObject *args, PyObject *kwds);

extern PyTypeObject ExtractCSVDataType;
