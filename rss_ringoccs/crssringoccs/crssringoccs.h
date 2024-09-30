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
#ifndef CRSSRINGOCCS_H
#define CRSSRINGOCCS_H

/*  The Python-C API is given here. The Python documentation recommends       *
 *  including Python.h before anything (even standard library headers).       */
#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif
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
typedef struct crssringoccs_PyDiffrecObj_Def {
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
} crssringoccs_PyDiffrecObj;

/*  The CSV struct containing all of the data for diffraction reconstruction. */
typedef struct crssringoccs_PyCSVObj_Def {
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
    PyObject *tau_phase_deg_vals;
    PyObject *tau_power_vals;
    PyObject *tau_vals;
} crssringoccs_PyCSVObj;

/******************************************************************************
 *  Function:                                                                 *
 *      crssringoccs_PyCSVObj_Destroy                                         *
 *  Purpose:                                                                  *
 *      Decrements all of the Python objects in the PyCSVObj struct. If no    *
 *      other references to the objects exists, they are free'd from memory.  *
 *  Arguments:                                                                *
 *      self (crssringoccs_PyCSVObj *):                                       *
 *          A pointer to the PyCSVObj that is to be destroyed.                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/py_csv_obj/                                 *
 *              crssringoccs_py_csv_obj_destroy.c                             *
 ******************************************************************************/
extern void crssringoccs_PyCSVObj_Destroy(crssringoccs_PyCSVObj *self);

/******************************************************************************
 *  Constant:                                                                 *
 *      crssringoccs_PyCSVObj_Members                                         *
 *  Purpose:                                                                  *
 *      Describes all of the atttributes in the PyCSVObj for the Python       *
 *      interpreter.                                                          *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/py_csv_obj/                                 *
 *              crssringoccs_py_csv_obj_members.c                             *
 ******************************************************************************/
extern PyMemberDef crssringoccs_PyCSVObj_Members[];

/******************************************************************************
 *  Constant:                                                                 *
 *      crssringoccs_PyCSVObj_Members                                         *
 *  Purpose:                                                                  *
 *      Describes all of the methods in the PyCSVObj for the Python           *
 *      interpreter.                                                          *
 *  Notes:                                                                    *
 *      PyCSVObj has no methods. This array only contains a NULL terminator.  *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/py_csv_obj/                                 *
 *              crssringoccs_py_csv_obj_members.c                             *
 ******************************************************************************/
extern PyMethodDef crssringoccs_PyCSVObj_Methods[];

/******************************************************************************
 *  Function:                                                                 *
 *      crssringoccs_ExtractCSVData_Steal                                     *
 *  Purpose:                                                                  *
 *      Steals data from a C struct and passes it to a Python object.         *
 *  Arguments:                                                                *
 *      py_csv (crssringoccs_PyCSVObj *):                                     *
 *          The Python object.                                                *
 *      csv (rssringoccs_CSVData *):                                          *
 *          The C struct containing the CSV data.                             *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/extract_csv_data/                           *
 *          crssringoccs_extract_csv_data_steal.c                             *
 ******************************************************************************/
extern void
crssringoccs_ExtractCSVData_Steal(crssringoccs_PyCSVObj *py_csv,
                                  rssringoccs_CSVData *csv);

/******************************************************************************
 *  Function:                                                                 *
 *      crssringoccs_ExtractCSVData_Init                                      *
 *  Purpose:                                                                  *
 *      Implements the __init__ method for the ExtractCSVData class.          *
 *  Arguments:                                                                *
 *      self (crssringoccs_PyCSVObj *):                                       *
 *          The Python object being initialized.                              *
 *      args (PyObject *):                                                    *
 *          The arguments to the ExtractCSVData class. These are the          *
 *          filenames to the CSV data.                                        *
 *      kwds (PyObject *):                                                    *
 *          The keywords to the ExtractCSVData class. These are the           *
 *          tau filename and the use_deprecate Boolean.                       *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/extract_csv_data/                           *
 *          crssringoccs_extract_csv_data_init.c                              *
 ******************************************************************************/
extern int
crssringoccs_ExtractCSVData_Init(crssringoccs_PyCSVObj *self,
                                 PyObject *args,
                                 PyObject *kwds);

/******************************************************************************
 *  Function:                                                                 *
 *      crssringoccs_ExtractCSVData_Create_History                            *
 *  Purpose:                                                                  *
 *      Creates the history dictionary for the ExtractCSVData class.          *
 *  Arguments:                                                                *
 *      self (crssringoccs_PyCSVObj *):                                       *
 *          The Python CSV object.                                            *
 *      geo_str (const char *):                                               *
 *          The path to the geo file.                                         *
 *      cal_str (const char *):                                               *
 *          The path to the cal file.                                         *
 *      dlp_str (const char *):                                               *
 *          The path to the dlp file.                                         *
 *      tau_str (const char *):                                               *
 *          The path to the tau file.                                         *
 *      use_deprecate (tmpl_Bool):                                            *
 *          The Boolean for using the older format.                           *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/extract_csv_data/                           *
 *          crssringoccs_extract_csv_data_create_history.c                    *
 ******************************************************************************/
extern void
crssringoccs_ExtractCSVData_Create_History(crssringoccs_PyCSVObj *self,
                                           const char *geo_str,
                                           const char *cal_str,
                                           const char *dlp_str,
                                           const char *tau_str,
                                           tmpl_Bool use_deprecate);

/******************************************************************************
 *  Constant:                                                                 *
 *      crssringoccs_ExtractCSVData                                           *
 *  Purpose:                                                                  *
 *      The ExtractCSVData class for python.                                  *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/extract_csv_data/                           *
 *              crssringoccs_get_uranus_data_class.c                          *
 ******************************************************************************/
extern PyTypeObject crssringoccs_ExtractCSVData;

/*  Data structure for the GEO.TAB files on the PDS.                          */
typedef struct PyGeoObj_Def {
    PyObject_HEAD
    PyObject *t_oet_spm_vals;
    PyObject *t_ret_spm_vals;
    PyObject *t_set_spm_vals;
    PyObject *rho_km_vals;
    PyObject *phi_rl_deg_vals;
    PyObject *phi_ora_deg_vals;
    PyObject *B_deg_vals;
    PyObject *D_km_vals;
    PyObject *rho_dot_kms_vals;
    PyObject *phi_rl_dot_kms_vals;
    PyObject *F_km_vals;
    PyObject *R_imp_km_vals;
    PyObject *rx_km_vals;
    PyObject *ry_km_vals;
    PyObject *rz_km_vals;
    PyObject *vx_kms_vals;
    PyObject *vy_kms_vals;
    PyObject *vz_kms_vals;
    PyObject *obs_spacecraft_lat_deg_vals;
    PyObject *history;
} rssringoccs_PyGeoObj;


extern void
crssringoccs_Create_Real_Numpy_Array(PyObject **py_ptr,
                                     double *ptr,
                                     size_t len);

extern void
crssringoccs_Create_Complex_Numpy_Array(PyObject **py_ptr,
                                        tmpl_ComplexDouble *ptr,
                                        size_t len);

extern void crssringoccs_Capsule_Cleanup(PyObject *capsule);

extern double *
crssringoccs_Extract_Data(rssringoccs_DLPObj *dlp,
                          PyObject *py_dlp,
                          const char *var_name);

extern rssringoccs_DLPObj *crssringoccs_Py_DLP_To_C_DLP(PyObject *py_dlp);

extern void
crssringoccs_C_Tau_To_Py_Tau(crssringoccs_PyDiffrecObj *py_tau,
                             rssringoccs_TAUObj *tau);

extern void
crssringoccs_Get_Py_Perturb(rssringoccs_TAUObj *tau, PyObject *perturb);

extern void
crssringoccs_Get_Py_Range(rssringoccs_TAUObj *tau, PyObject *rngreq);

extern void
crssringoccs_Get_Py_Vars_From_Tau_Self(rssringoccs_TAUObj *tau,
                                       const crssringoccs_PyDiffrecObj *self);

extern PyMemberDef crssringoccs_DiffractionCorrection_Members[];

extern PyMethodDef crssringoccs_DiffractionCorrection_Methods[];

extern void
crssringoccs_DiffractionCorrection_Destroy(crssringoccs_PyDiffrecObj *self);

extern int
crssringoccs_DiffractionCorrection_Init(crssringoccs_PyDiffrecObj *self,
                                        PyObject *args,
                                        PyObject *kwds);

extern PyTypeObject crssringoccs_DiffractionCorrection;

extern void
crssringoccs_GetUranusData_Steal(crssringoccs_PyCSVObj *py_csv,
                                 rssringoccs_UranusCSVData *csv);

extern int
crssringoccs_GetUranusData_Init(crssringoccs_PyCSVObj *self,
                                PyObject *args,
                                PyObject *kwds);

extern void
crssringoccs_GetUranusData_Create_History(crssringoccs_PyCSVObj *self,
                                          const char *geo_str,
                                          const char *dlp_str,
                                          const char *tau_str,
                                          tmpl_Bool dlp_in_radians);

extern PyTypeObject crssringoccs_GetUranusData;

extern void
crssringoccs_GetMergedCSVData_Steal(crssringoccs_PyCSVObj *py_csv,
                                    rssringoccs_MergedCSVData *csv);

extern int
crssringoccs_GetMergedCSVData_Init(crssringoccs_PyCSVObj *self,
                                   PyObject *args,
                                   PyObject *kwds);

extern void
crssringoccs_GetMergedCSVData_Create_History(crssringoccs_PyCSVObj *self,
                                             const char *dlpm_str);

extern PyTypeObject crssringoccs_GetMergedCSVData;

#endif
