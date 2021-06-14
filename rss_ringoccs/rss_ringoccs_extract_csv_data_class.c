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
 ******************************************************************************
 *                        Diffraction Correction Class                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Defines the DiffractionCorrection class for rss_ringoccs. This uses   *
 *      the C-Python API to build an extension module containing the class.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 22, 2019                                                 *
 ******************************************************************************/

/*  To avoid compiler warnings about deprecated numpy stuff.                  */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN

/*  The following are NON-STANDARD header files for using the C-Python API,   *
 *  and the numpy API.                                                        */
#include <Python.h>
#include <structmember.h>
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

/*  The standard library header stdlib contains malloc, calloc, and realloc.  */
#include <stdlib.h>

/*  The puts function is found here.                                          */
#include <stdio.h>

/*  The following header files come from libtmpl. This is a NON-STANDARD      *
 *  library that contains a variety of math tools.                            */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>

/*  And a bunch of headers from this project.                                 */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  The definition of the DiffractionCorrection class as a C struct.          */
typedef struct _PyCSVObj {
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

/*  This function frees the memory allocated to a pointer by malloc when the  *
 *  corresponding variable is destroyed at the Python level. Without this you *
 *  will have serious memory leaks, so do not remove!                         */
static void capsule_cleanup(PyObject *capsule)
{
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    free(memory);
}

static void
set_var(PyObject **py_ptr, double **ptr, unsigned long int len)
{
    PyObject *arr, *tmp, *capsule;
    long pylength = (long)len;

    arr = PyArray_SimpleNewFromData(1, &pylength, NPY_DOUBLE, *ptr);
    capsule = PyCapsule_New((void *) (*ptr), NULL, capsule_cleanup);

    PyArray_SetBaseObject((PyArrayObject *)arr, capsule);

    tmp = *py_ptr;
    Py_INCREF(arr);
    *py_ptr = arr;
    Py_XDECREF(tmp);
}

static void
rssringoccs_C_CSV_to_Py_CSV(PyCSVObj *py_csv, rssringoccs_CSVData *csv)
{
    if (csv == NULL)
        return;

    if (csv->error_occurred)
        return;

    if (py_csv == NULL)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_CSV_to_Py_CSV\n\n"
            "\rInput py_csv is NULL. Aborting.n"
        );
        return;
    }

    set_var(&py_csv->rho_km_vals, &csv->rho_km_vals, csv->n_elements);
    set_var(&py_csv->B_rad_vals, &csv->B_rad_vals, csv->n_elements);
    set_var(&py_csv->D_km_vals, &csv->D_km_vals, csv->n_elements);
    set_var(&py_csv->f_sky_hz_vals, &csv->f_sky_hz_vals, csv->n_elements);
    set_var(&py_csv->p_norm_vals, &csv->p_norm_vals, csv->n_elements);
    set_var(&py_csv->phase_rad_vals, &csv->phase_rad_vals, csv->n_elements);
    set_var(&py_csv->phi_rad_vals, &csv->phi_rad_vals, csv->n_elements);
    set_var(&py_csv->phi_rl_rad_vals, &csv->phi_rl_rad_vals, csv->n_elements);
    set_var(&py_csv->rho_dot_kms_vals, &csv->rho_dot_kms_vals, csv->n_elements);
    set_var(&py_csv->t_oet_spm_vals, &csv->t_oet_spm_vals, csv->n_elements);
    set_var(&py_csv->t_ret_spm_vals, &csv->t_ret_spm_vals, csv->n_elements);
    set_var(&py_csv->t_set_spm_vals, &csv->t_set_spm_vals, csv->n_elements);
    set_var(&py_csv->tau_vals, &csv->tau_vals, csv->n_elements);
    set_var(&py_csv->rx_km_vals, &csv->rx_km_vals, csv->n_elements);
    set_var(&py_csv->ry_km_vals, &csv->ry_km_vals, csv->n_elements);
    set_var(&py_csv->rz_km_vals, &csv->rz_km_vals, csv->n_elements);
    set_var(&py_csv->raw_tau_threshold_vals, &csv->raw_tau_threshold_vals, csv->n_elements);
    set_var(&py_csv->rho_corr_pole_km_vals, &csv->rho_corr_pole_km_vals, csv->n_elements);
    set_var(&py_csv->rho_corr_timing_km_vals, &csv->rho_corr_timing_km_vals, csv->n_elements);
}

/*  Deallocating function for the DiffractionCorrection class.                */
static void ExtractCSVData_dealloc(PyCSVObj *self)
{
    Py_XDECREF(self->B_rad_vals);
    Py_XDECREF(self->D_km_vals);
    Py_XDECREF(self->f_sky_hz_vals);
    Py_XDECREF(self->p_norm_vals);
    Py_XDECREF(self->phase_rad_vals);
    Py_XDECREF(self->phi_rad_vals);
    Py_XDECREF(self->phi_rl_rad_vals);
    Py_XDECREF(self->raw_tau_threshold_vals);
    Py_XDECREF(self->rev_info);
    Py_XDECREF(self->rho_corr_pole_km_vals);
    Py_XDECREF(self->rho_corr_timing_km_vals);
    Py_XDECREF(self->rho_dot_kms_vals);
    Py_XDECREF(self->rho_km_vals);
    Py_XDECREF(self->t_oet_spm_vals);
    Py_XDECREF(self->t_ret_spm_vals);
    Py_XDECREF(self->t_set_spm_vals);
    Py_XDECREF(self->tau_vals);
    Py_XDECREF(self->history);
    Py_XDECREF(self->input_vars);
    Py_XDECREF(self->input_kwds);
    Py_XDECREF(self->rx_km_vals);
    Py_XDECREF(self->ry_km_vals);
    Py_XDECREF(self->rz_km_vals);
    Py_XDECREF(self->tau_rho);
    Py_XDECREF(self->tau_phase);
    Py_XDECREF(self->tau_power);
    Py_XDECREF(self->tau_vals);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

/*  The init function for the dirrection correction class. This is the        *
 *  equivalent of the __init__ function defined in a normal python class.     */
static int ExtractCSVData_init(PyCSVObj *self, PyObject *args, PyObject *kwds)
{
    rssringoccs_CSVData *csv;
    tmpl_Bool dpr = tmpl_False;

    /*  The list of the keywords accepted by the DiffractionCorrection class. *
     *  dlp and res are REQUIRED inputs, the rest are optional. If the user   *
     *  does not provide these optional keywords, we must set them ourselves. */
    static char *kwlist[] = {"geo", "cal", "dlp", "use_deprecate", "tau", NULL};

    /*  Python objects needed throughout the computation.                     */
    PyObject *tmp, *csv_tmp;
    const char *geo_str = NULL;
    const char *cal_str = NULL;
    const char *dlp_str = NULL;
    const char *tau_str = NULL;

    /*  Extract the inputs and keywords supplied by the user. If the data     *
     *  cannot be extracted, raise a type error and return to caller. A short *
     *  explaination of PyArg_ParseTupleAndKeywords. The inputs args and kwds *
     *  are somewhat straight-forward, they're the arguments and keywords     *
     *  passed by the string. The cryptic string is not straight-forward. The *
     *  | symbol means everything after need not be positional, and we can    *
     *  specify arguments and keywords by name when calling                   *
     *  DiffractionCorrection, for example                                    *
     *  DiffractionCorrect(..., wtype="blah"). O indicates a Python object,   *
     *  and d is a Python float. This is the DLP and res variables. The $     *
     *  symbold means everything after is optional. s is a string, p is a     *
     *  Boolean (p for "predicate"). b is an integer, and the colon : denotes *
     *  that the input list has ended.                                        */
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|sss$ps:", kwlist, &geo_str,
                                     &cal_str, &dlp_str, &dpr, &tau_str))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcsvtools.ExtractCSVData\n\n"
            "\rCould not parse input variables.\n\n"
            "\rInputs:\n"
            "\r\tgeo:           Location of a GEO.TAB file (str)\n"
            "\r\tcal:           Location of a CAL.TAB file (str)\n"
            "\r\tdlp:           Location of a DLP.TAB file (str)\n"
            "\rKeywords:\n"
            "\r\tuse_deprecate: Use the old CSV format.\n"
            "\r\ttau:           Location of a TAU.TAB file (str)\n"
        );
        return -1;
    }

    csv = rssringoccs_Extract_CSV_Data(geo_str, cal_str, dlp_str, tau_str, dpr);

    if (csv == NULL)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rFailed to pass variables to C. rssringoccs_Py_DLP_to_C_DLP\n"
            "\rreturned NULL. Returning.\n\n"
        );
        rssringoccs_Destroy_CSV(&csv);
        return -1;
    }

    if (csv->error_occurred)
    {
        if (csv->error_message == NULL)
        {
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rFailed to pass variables to C. rssringoccs_Py_DLP_to_C_DLP\n"
                "\rreturned a dlp with error_occurred set to True. No\n"
                "\rerror message was set. Returning.\n\n"
            );
        }
        else
            PyErr_Format(PyExc_RuntimeError, "%s", csv->error_message);

        rssringoccs_Destroy_CSV(&csv);
        return -1;
    }

    rssringoccs_C_CSV_to_Py_CSV(self, csv);

    csv_tmp = Py_BuildValue(
        "{s:s,s:s,s:s}",
        "geo", geo_str,
        "cal", cal_str,
        "dlp", dlp_str
    );

    tmp = self->input_vars;
    Py_INCREF(csv_tmp);
    self->input_vars = csv_tmp;
    Py_XDECREF(tmp);

    return 1;
}

static PyMemberDef Custom_members[] = {
    {
        "rho_km_vals", T_OBJECT_EX, offsetof(PyCSVObj, rho_km_vals), 0,
        "Ring radius."
    },
    {
        "phase_rad_vals", T_OBJECT_EX, offsetof(PyCSVObj, phase_rad_vals),
        0, "Raw diffracted phase."
    },
    {
        "B_rad_vals", T_OBJECT_EX, offsetof(PyCSVObj, B_rad_vals), 0,
        "Ring inclination angle."
    },
    {
        "D_km_vals", T_OBJECT_EX, offsetof(PyCSVObj, D_km_vals), 0,
        "Spacecraft to ring-intercept point distance."
    },
    {
        "rx_km_vals", T_OBJECT_EX, offsetof(PyCSVObj, rx_km_vals), 0,
        "x coordinate of the spacecraft in planetocentric frame."
    },
    {
        "ry_km_vals", T_OBJECT_EX, offsetof(PyCSVObj, ry_km_vals), 0,
        "y coordinate of the spacecraft in planetocentric frame."
    },
    {
        "rz_km_vals", T_OBJECT_EX, offsetof(PyCSVObj, rz_km_vals), 0,
        "z coordinate of the spacecraft in planetocentric frame."
    },
    {
        "f_sky_hz_vals", T_OBJECT_EX, offsetof(PyCSVObj, f_sky_hz_vals), 0,
        "Frequency of the input signal"
    },
    {
        "p_norm_vals", T_OBJECT_EX, offsetof(PyCSVObj, p_norm_vals), 0,
        "Raw power data"
    },
    {
        "phi_rad_vals", T_OBJECT_EX, offsetof(PyCSVObj, phi_rad_vals), 0,
        "Ring azimuth angle"
    },
    {
        "phi_rl_rad_vals", T_OBJECT_EX, offsetof(PyCSVObj, phi_rl_rad_vals),
        0, "Ring longitude angle"
    },
    {
        "raw_tau_threshold_vals", T_OBJECT_EX,
        offsetof(PyCSVObj, raw_tau_threshold_vals), 0,
        "Raw threshold optical depth"
    },
    {
        "rev_info", T_OBJECT_EX, offsetof(PyCSVObj, rev_info), 0,
        "Information about the occultation"
    },
    {
        "rho_corr_pole_km_vals", T_OBJECT_EX,
        offsetof(PyCSVObj, rho_corr_pole_km_vals), 0,
        "Ring radius with pole correction."
    },
    {
        "rho_corr_timing_km_vals", T_OBJECT_EX,
        offsetof(PyCSVObj, rho_corr_timing_km_vals), 0,
        "Ring radius with timing correction."
    },
    {
        "rho_dot_kms_vals", T_OBJECT_EX,
        offsetof(PyCSVObj, rho_dot_kms_vals), 0,
        "Time derivative of the ring radius."
    },
    {
        "t_oet_spm_vals", T_OBJECT_EX, offsetof(PyCSVObj, t_oet_spm_vals),
        0, "Observed event time in seconds past midnight"
    },
    {
        "t_ret_spm_vals", T_OBJECT_EX, offsetof(PyCSVObj, t_ret_spm_vals),
        0, "Ring event time in seconds past midnight"
    },
    {
        "t_set_spm_vals", T_OBJECT_EX, offsetof(PyCSVObj, t_set_spm_vals),
        0, "Spacecraft event time in seconds past midnight"
    },
    {
        "tau_vals", T_OBJECT_EX, offsetof(PyCSVObj, tau_vals), 0,
        "Optical depth"
    },
    {
        "history", T_OBJECT_EX, offsetof(PyCSVObj, history), 0,
        "History of the tau instance"
    },
    {
        "input_vars", T_OBJECT_EX, offsetof(PyCSVObj, input_vars), 0,
        "Dictionary of input arguments used to create this instance."
    },
    {
        "input_kwds", T_OBJECT_EX, offsetof(PyCSVObj, input_kwds), 0,
        "Dictionary of input keywords used to create this instance."
    },
    {
        NULL
    }  /* Sentinel */
};

static PyObject *
ExtractCSVData(PyCSVObj *self, PyObject *Py_UNUSED(ignored))
{
    return PyUnicode_FromFormat("ExtractCSVData");
}

static PyMethodDef ExtractCSVData_methods[] =
{
    {
        "ExtractCSVData",
        (PyCFunction) ExtractCSVData,
        METH_NOARGS,
        "ExtractCSVData class."
    },
    {NULL}
};

static PyTypeObject ExtractCSVDataType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "csvtools.ExtractCSVData",
    .tp_doc = "ExtractCSVData class.",
    .tp_basicsize = sizeof(PyCSVObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) ExtractCSVData_init,
    .tp_dealloc = (destructor) ExtractCSVData_dealloc,
    .tp_members = Custom_members,
    .tp_methods = ExtractCSVData_methods,
};

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "custom",
    .m_doc = "Module containing the rss_ringoccs class.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_csvtools(void)
{
    PyObject *m;
    int pymod_addobj;
    if (PyType_Ready(&ExtractCSVDataType) < 0)
        return NULL;

    m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    Py_INCREF(&ExtractCSVDataType);
    pymod_addobj = PyModule_AddObject(m, "ExtractCSVData",
                                      (PyObject *) &ExtractCSVDataType);

    if (pymod_addobj < 0)
    {
        Py_DECREF(&ExtractCSVDataType);
        Py_DECREF(m);
        return NULL;
    }

    import_array();
    return m;
}
