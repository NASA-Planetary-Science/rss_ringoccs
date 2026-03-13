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
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Macro for setting a variable to the Python "None" type.                   */
#define MAKE_NONE(var)                                                         \
    do {                                                                       \
        PyObject *tmp = var;                                                   \
        Py_INCREF(Py_None);                                                    \
        var = Py_None;                                                         \
        Py_CLEAR(tmp);                                                         \
    } while (0)

/*  The definition of the DiffractionCorrection class as a C struct.          */
typedef struct crssringoccs_PyDiffrecObj_Def {
    PyObject_HEAD
    rssringoccs_TAUObj *tau;          /*  Tau object with all of the C data.  */
    PyObject *T_in;                   /*  Input complex transmittance.        */
    PyObject *T_out;                  /*  Reconstructed complex transmittance.*/
    PyObject *T_fwd;                  /*  Forward model complex transmittance.*/
    PyObject *k_vals;                 /*  Wavenumber, 2 pi / wavelength.      */
    PyObject *B_deg_vals;             /*  Ring opening angle.                 */
    PyObject *D_km_vals;              /*  Spacecraft-to-Ring distance.        */
    PyObject *F_km_vals;              /*  Fresnel scale.                      */
    PyObject *phi_deg_vals;           /*  Ring azimuth angle.                 */
    PyObject *phi_rl_deg_vals;        /*  Ring longitude angle.               */
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
    PyObject *input_vars;             /*  Input parameters for the class.     */
    PyObject *input_kwds;             /*  Input keywords for the class.       */
    PyObject *rngreq;                 /*  Requested range, Python keyword.    */
    PyObject *perturb;                /*  Polynomial perturbation for psi.    */

    /*  The following members are created with setter functions at the Python *
     *  level. The are initialized to NULL and left alone unless the user     *
     *  explicitly requests them. All of these are computable from the other  *
     *  members in this struct, hence we avoid computing them unless this is  *
     *  explicitly desired to avoid wasting memory.                           */
    PyObject *p_norm_vals;            /*  Normalized DLP power values.        */
    PyObject *power_vals;             /*  Reconstructed power values.         */
    PyObject *p_fwd_vals;             /*  Forward model power values.         */
    PyObject *phase_norm_deg_vals;    /*  Normalized DLP phase values.        */
    PyObject *phase_deg_vals;         /*  Reconstructed phase values.         */
    PyObject *phase_fwd_deg_vals;     /*  Forward model phase values.         */
    PyObject *tau_norm_vals;          /*  Normalized DLP optical depth values.*/
    PyObject *tau_vals;               /*  Reconstructed optical depth values. */
    PyObject *tau_fwd_vals;           /*  Forward model optical depth values. */

    /*  The remaining members are not arrays, but hold the keywords and       *
     *  arguments passed to the DiffractionCorrection class when initialized. */
    tmpl_Bool bfac;                   /*  Boolean for b factor in resolution. */
    tmpl_Bool use_fwd;                /*  Boolean for forward modeling.       */
    tmpl_Bool use_norm;               /*  Boolean for window normalization.   */
    tmpl_Bool verbose;                /*  Boolean for printing messages.      */
    double input_resolution_km;       /*  Input resolution, in kilometers.    */
    double resolution_factor;         /*  Resolution scale factor, unitless.  */
    double eccentricity;              /*  Eccentricity, elliptical rings only.*/
    double periapse;                  /*  Periapse, elliptical rings only.    */
    double sigma;                     /*  Allen deviation of spacecraft.      */
    const char *outfiles;             /*  TAB files for this Tau object.      */
    const char *wtype;                /*  Requested window type.              */
    const char *psitype;              /*  Requested reconstruction algorithm. */
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
 *      crssringoccs_CassiniCSVData_Steal                                     *
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
 *      rss_ringoccs/crssringoccs/cassini_csv_data/                           *
 *          crssringoccs_cassini_csv_data_steal.c                             *
 ******************************************************************************/
extern void
crssringoccs_CassiniCSVData_Steal(crssringoccs_PyCSVObj *py_csv,
                                  rssringoccs_CSVData *csv);

/******************************************************************************
 *  Function:                                                                 *
 *      crssringoccs_CassiniCSVData_Init                                      *
 *  Purpose:                                                                  *
 *      Implements the __init__ method for the CassiniCSVData class.          *
 *  Arguments:                                                                *
 *      self (crssringoccs_PyCSVObj *):                                       *
 *          The Python object being initialized.                              *
 *      args (PyObject *):                                                    *
 *          The arguments to the CassiniCSVData class. These are the          *
 *          filenames to the CSV data.                                        *
 *      kwds (PyObject *):                                                    *
 *          The keywords to the CassiniCSVData class. These are the           *
 *          tau filename and the use_deprecate Boolean.                       *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/cassini_csv_data/                           *
 *          crssringoccs_cassini_csv_data_init.c                              *
 ******************************************************************************/
extern int
crssringoccs_CassiniCSVData_Init(crssringoccs_PyCSVObj *self,
                                 PyObject *args,
                                 PyObject *kwds);

/******************************************************************************
 *  Function:                                                                 *
 *      crssringoccs_CassiniCSVData_Create_History                            *
 *  Purpose:                                                                  *
 *      Creates the history dictionary for the CassiniCSVData class.          *
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
 *      rss_ringoccs/crssringoccs/cassini_csv_data/                           *
 *          crssringoccs_cassini_csv_data_create_history.c                    *
 ******************************************************************************/
extern void
crssringoccs_CassiniCSVData_Create_History(crssringoccs_PyCSVObj *self,
                                           const char *geo_str,
                                           const char *cal_str,
                                           const char *dlp_str,
                                           const char *tau_str,
                                           tmpl_Bool use_deprecate);

/******************************************************************************
 *  Constant:                                                                 *
 *      crssringoccs_CassiniCSVData                                           *
 *  Purpose:                                                                  *
 *      The CassiniCSVData class for python.                                  *
 *  Source Code:                                                              *
 *      rss_ringoccs/crssringoccs/cassini_csv_data/                           *
 *              crssringoccs_cassini_csv_data_class.c                         *
 ******************************************************************************/
extern PyTypeObject crssringoccs_CassiniCSVData;

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
crssringoccs_Create_Real_Numpy_Array(PyObject ** const py_ptr,
                                     double * const ptr,
                                     void (*cleanup)(PyObject *),
                                     const size_t len);

extern void
crssringoccs_Create_Complex_Numpy_Array(PyObject ** const py_ptr,
                                        tmpl_ComplexDouble * const ptr,
                                        void (*cleanup)(PyObject *),
                                        const size_t len);

extern void crssringoccs_Capsule_Cleanup(PyObject * const capsule);

extern double *
crssringoccs_DLP_Extract_Data(rssringoccs_DLPObj * const dlp,
                              PyObject * const object,
                              const char * const var_name);

extern rssringoccs_DLPObj *
crssringoccs_PyObject_To_DLP(PyObject * const object);

extern void
crssringoccs_DiffractionCorrection_Set_Perturb(
    crssringoccs_PyDiffrecObj * const self
);

extern void
crssringoccs_DiffractionCorrection_Set_Range(
    crssringoccs_PyDiffrecObj * const self
);

extern void
crssringoccs_DiffractionCorrection_Set_Keywords(
    crssringoccs_PyDiffrecObj * const self
);

extern void
crssringoccs_DiffractionCorrection_Create_Argument_Dictionary(
    crssringoccs_PyDiffrecObj * const self,
    PyObject * const dlp
);

extern void
crssringoccs_DiffractionCorrection_Create_Keyword_Dictionary(
    crssringoccs_PyDiffrecObj * const self
);

extern void
crssringoccs_DiffractionCorrection_Finish(
    crssringoccs_PyDiffrecObj * const self, PyObject * const dlp
);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_P_Norm_Vals(PyObject *op, void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_Power_Vals(PyObject *op, void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_P_Fwd_Vals(PyObject *op, void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_P_Norm_Vals(PyObject *op,
                                                   PyObject *value,
                                                   void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_Power_Vals(PyObject *op,
                                                  PyObject *value,
                                                  void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_P_Fwd_Vals(PyObject *op,
                                                  PyObject *value,
                                                  void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_Phase_Norm_Deg_Vals(PyObject *op,
                                                           void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_Phase_Deg_Vals(PyObject *op,
                                                      void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_Phase_Fwd_Deg_Vals(PyObject *op,
                                                          void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_Phase_Norm_Deg_Vals(PyObject *op,
                                                           PyObject *value,
                                                           void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_Phase_Deg_Vals(PyObject *op,
                                                      PyObject *value,
                                                      void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_Phase_Fwd_Deg_Vals(PyObject *op,
                                                          PyObject *value,
                                                          void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_Tau_Vals(PyObject *op,
                                                void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_Tau_Norm_Vals(PyObject *op,
                                                     void *closure);

extern PyObject *
crssringoccs_DiffractionCorrection_Get_Tau_Fwd_Vals(PyObject *op,
                                                    void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_Tau_Vals(PyObject *op,
                                                PyObject *value,
                                                void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_Tau_Norm_Vals(PyObject *op,
                                                     PyObject *value,
                                                     void *closure);

extern int
crssringoccs_DiffractionCorrection_Set_Tau_Fwd_Vals(PyObject *op,
                                                    PyObject *value,
                                                    void *closure);

extern void
crssringoccs_DiffractionCorrection_Destroy(crssringoccs_PyDiffrecObj *self);

extern PyObject *
crssringoccs_DiffractionCorrection_New(PyTypeObject *type,
                                       PyObject *args,
                                       PyObject *kwds);

extern int
crssringoccs_DiffractionCorrection_Init(crssringoccs_PyDiffrecObj *self,
                                        PyObject *args,
                                        PyObject *kwds);

extern PyGetSetDef crssringoccs_DiffractionCorrection_GetSetters[];
extern PyMemberDef crssringoccs_DiffractionCorrection_Members[];
extern PyMethodDef crssringoccs_DiffractionCorrection_Methods[];
extern PyTypeObject crssringoccs_DiffractionCorrection;

extern void
crssringoccs_VoyagerCSVData_Steal(crssringoccs_PyCSVObj *py_csv,
                                  rssringoccs_UranusCSVData *csv);

extern int
crssringoccs_VoyagerCSVData_Init(crssringoccs_PyCSVObj *self,
                                 PyObject *args,
                                 PyObject *kwds);

extern void
crssringoccs_VoyagerCSVData_Create_History(crssringoccs_PyCSVObj *self,
                                           const char *geo_str,
                                           const char *dlp_str,
                                           const char *tau_str,
                                           tmpl_Bool dlp_in_radians);

extern PyTypeObject crssringoccs_VoyagerCSVData;

extern void
crssringoccs_MergedCSVData_Steal(crssringoccs_PyCSVObj *py_csv,
                                 rssringoccs_MergedCSVData *csv);

extern int
crssringoccs_MergedCSVData_Init(crssringoccs_PyCSVObj *self,
                                PyObject *args,
                                PyObject *kwds);

extern void
crssringoccs_MergedCSVData_Create_History(crssringoccs_PyCSVObj *self,
                                          const char *dlpm_str);

extern PyTypeObject crssringoccs_MergedCSVData;

#endif
