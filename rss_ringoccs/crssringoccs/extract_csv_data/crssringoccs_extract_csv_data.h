#ifndef CRSSRINGOCCS_EXTRACT_CSV_DATA_H
#define CRSSRINGOCCS_EXTRACT_CSV_DATA_H

/*  The Python-C API is given here. The Python documentation recommends       *
 *  including Python.h before anything (even standard library headers).       */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

/*  And a bunch of headers from this project.                                 */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

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
} crssringoccs_PyCSVObj;

extern void
crssringoccs_ExtractCSVData_Steal(crssringoccs_PyCSVObj *py_csv,
                                  rssringoccs_CSVData *csv);

extern void crssringoccs_ExtractCSVData_Destroy(crssringoccs_PyCSVObj *self);

extern int
crssringoccs_ExtractCSVData_Init(crssringoccs_PyCSVObj *self,
                                 PyObject *args,
                                 PyObject *kwds);

extern PyTypeObject ExtractCSVDataType;

#endif
