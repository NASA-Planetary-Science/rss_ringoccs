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
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/tmpl_special_functions.h>

/*  And a bunch of headers from this project.                                 */
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  The definition of the DiffractionCorrection class as a C struct.          */
typedef struct _PyDiffrecObj {
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

static double *
extract_data(rssringoccs_DLPObj *dlp, PyObject *py_dlp, const char *var_name)
{
    PyObject *tmp;
    PyObject *arr;
    unsigned long len;

    if (dlp == NULL)
        return NULL;

    if (dlp->error_occurred)
        return NULL;

    if (py_dlp == NULL)
        return NULL;

    if (!PyObject_HasAttrString(py_dlp, var_name))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\t%s\n\n",
            var_name
        );
        return NULL;
    }
    else
        tmp = PyObject_GetAttrString(py_dlp, var_name);

    if (!PyArray_Check(tmp))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a numpy array.\n",
            var_name
        );
        return NULL;
    }
    else
        arr = PyArray_FromObject(tmp, NPY_DOUBLE, 1, 1);

    len = (unsigned long)PyArray_DIMS((PyArrayObject *)arr)[0];

    /*  If PyArray_FromObject failed arr should be NULL. If so, raise error.  */
    if (!arr)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a numpy array.\n",
            var_name
        );
        return NULL;
    }

    /*  Currently we only allow for one dimensional inputs.                   */
    else if (PyArray_NDIM((PyArrayObject *)arr) != 1)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s must be a one-dimensional numpy array.\n",
            var_name
        );
        return NULL;
    }

    /*  arr should have the same number of elements as rho_km_vals.           */
    else if (len != dlp->arr_size)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = malloc(sizeof(*dlp->error_message) * 256);
        sprintf(
            dlp->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\r%s and rho_km_vals have a different number of elements.\n",
            var_name
        );
        return NULL;
    }

    /*  If every passed, set ptr to point to the data inside the array arr.   */
    return (double *)PyArray_DATA((PyArrayObject *)arr);
}

static rssringoccs_DLPObj *
rssringoccs_Py_DLP_to_C_DLP(PyObject *py_dlp)
{
    PyObject *tmp;
    PyObject *arr;
    rssringoccs_DLPObj *dlp;

    if (py_dlp == NULL)
        return NULL;

    dlp = malloc(sizeof(*dlp));
    if (dlp == NULL)
        return dlp;

    dlp->error_occurred = tmpl_False;
    dlp->error_message = NULL;

    if (py_dlp == NULL)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is NULL.\n"
        );
        return dlp;
    }

    /*  Next we're going to run error checks on the input numpy arrays which  *
     *  should be contained inside of the DLPInst object. We'll check that    *
     *  these attributes exist, that they are numpy arrays, are 1 dimensional,*
     *  and have the same number of elements as rho_km_vals. We'll also       *
     *  convert the arrays to double and retrieve a pointer to the data.      *
     *  First, we need to make sure rho_km_vals is a legal numpy array and    *
     *  extract the length of it. Check that rho_km_vals exists in DLPInst.   */
    if (!PyObject_HasAttrString(py_dlp, "rho_km_vals"))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\trho_km_vals\n\n"
        );
        return dlp;
    }

    /*  If it exists, get a pointer to it.                                    */
    else
        tmp = PyObject_GetAttrString(py_dlp, "rho_km_vals");

    /*  Now make sure rho_km_vals is a numpy array.                           */
    if (!PyArray_Check(tmp))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrho_km_vals must be a numpy array.\n"
        );
        return dlp;
    }

    /*  If rho_km_vals is a numpy array, try to convert it to double.         */
    else
        arr = PyArray_FromObject(tmp, NPY_DOUBLE, 1, 1);

    /*  If PyArray_FromObject failed arr should be NULL. If so, raise error. */
    if (!arr)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rCould not convert rho_km_vals to double array. Input is most\n"
            "\rlikely complex numbers or contains a string.\n\n"
        );
        return dlp;
    }

    /*  Currently we only allow for one dimensional inputs.                   */
    else if (PyArray_NDIM((PyArrayObject *)arr) != 1)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrho_km_vals must be a one-dimensional numpy array.\n"
        );
        return dlp;
    }

    /*  If every passed, set tau.rho_km_vals to point to the data inside arr. */
    dlp->rho_km_vals = (double *)PyArray_DATA((PyArrayObject *)arr);
    dlp->arr_size = PyArray_DIMS((PyArrayObject *)arr)[0];

    dlp->p_norm_vals = extract_data(dlp, py_dlp, "p_norm_vals");
    dlp->phase_rad_vals = extract_data(dlp, py_dlp, "phase_rad_vals");
    dlp->phi_rad_vals = extract_data(dlp, py_dlp, "phi_rad_vals");
    dlp->phi_rl_rad_vals = extract_data(dlp, py_dlp, "phi_rl_rad_vals");
    dlp->B_rad_vals = extract_data(dlp, py_dlp, "B_rad_vals");
    dlp->D_km_vals = extract_data(dlp, py_dlp, "D_km_vals");
    dlp->f_sky_hz_vals = extract_data(dlp, py_dlp, "f_sky_hz_vals");
    dlp->rho_dot_kms_vals = extract_data(dlp, py_dlp, "rho_dot_kms_vals");
    dlp->t_oet_spm_vals = extract_data(dlp, py_dlp, "t_oet_spm_vals");
    dlp->t_ret_spm_vals = extract_data(dlp, py_dlp, "t_ret_spm_vals");
    dlp->t_set_spm_vals = extract_data(dlp, py_dlp, "t_set_spm_vals");
    dlp->rx_km_vals = extract_data(dlp, py_dlp, "rx_km_vals");
    dlp->ry_km_vals = extract_data(dlp, py_dlp, "ry_km_vals");
    dlp->rz_km_vals = extract_data(dlp, py_dlp, "rz_km_vals");
    dlp->rho_corr_pole_km_vals = extract_data(dlp, py_dlp, "rho_corr_pole_km_vals");
    dlp->rho_corr_timing_km_vals = extract_data(dlp, py_dlp, "rho_corr_timing_km_vals");
    dlp->raw_tau_threshold_vals = extract_data(dlp, py_dlp, "raw_tau_threshold_vals");
    return dlp;
}

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
    PyObject *arr;
    PyObject *capsule;
    PyObject *tmp;
    long pylength = (long)len;

    arr     = PyArray_SimpleNewFromData(1, &pylength, NPY_DOUBLE, *ptr);
    capsule = PyCapsule_New((void *) (*ptr), NULL, capsule_cleanup);

    PyArray_SetBaseObject((PyArrayObject *)arr, capsule);

    tmp = *py_ptr;
    Py_INCREF(arr);
    *py_ptr = arr;
    Py_XDECREF(tmp);
}

static void rssringoccs_C_Tau_to_Py_Tau(PyDiffrecObj *py_tau,
                                        rssringoccs_TAUObj *tau)
{
    PyObject *tmp;
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (py_tau == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_Tau_to_Py_Tau\n\n"
            "\rInput py_tau is NULL. Aborting.n"
        );
        return;
    }

    set_var(&py_tau->rho_km_vals, &tau->rho_km_vals, tau->arr_size);
    set_var(&py_tau->B_rad_vals, &tau->B_rad_vals, tau->arr_size);
    set_var(&py_tau->D_km_vals, &tau->D_km_vals, tau->arr_size);
    set_var(&py_tau->F_km_vals, &tau->F_km_vals, tau->arr_size);
    set_var(&py_tau->f_sky_hz_vals, &tau->f_sky_hz_vals, tau->arr_size);
    set_var(&py_tau->p_norm_vals, &tau->p_norm_vals, tau->arr_size);
    set_var(&py_tau->phase_rad_vals, &tau->phase_rad_vals, tau->arr_size);
    set_var(&py_tau->phase_vals, &tau->phase_vals, tau->arr_size);
    set_var(&py_tau->phi_rad_vals, &tau->phi_rad_vals, tau->arr_size);
    set_var(&py_tau->phi_rl_rad_vals, &tau->phi_rl_rad_vals, tau->arr_size);
    set_var(&py_tau->power_vals, &tau->power_vals, tau->arr_size);
    set_var(&py_tau->rho_dot_kms_vals, &tau->rho_dot_kms_vals, tau->arr_size);
    set_var(&py_tau->t_oet_spm_vals, &tau->t_oet_spm_vals, tau->arr_size);
    set_var(&py_tau->t_ret_spm_vals, &tau->t_ret_spm_vals, tau->arr_size);
    set_var(&py_tau->t_set_spm_vals, &tau->t_set_spm_vals, tau->arr_size);
    set_var(&py_tau->tau_vals, &tau->tau_vals, tau->arr_size);
    set_var(&py_tau->w_km_vals, &tau->w_km_vals, tau->arr_size);
    set_var(&py_tau->rx_km_vals, &tau->rx_km_vals, tau->arr_size);
    set_var(&py_tau->ry_km_vals, &tau->ry_km_vals, tau->arr_size);
    set_var(&py_tau->rz_km_vals, &tau->rz_km_vals, tau->arr_size);
    set_var(&py_tau->raw_tau_threshold_vals, &tau->raw_tau_threshold_vals, tau->arr_size);
    set_var(&py_tau->rho_corr_pole_km_vals, &tau->rho_corr_pole_km_vals, tau->arr_size);
    set_var(&py_tau->rho_corr_timing_km_vals, &tau->rho_corr_timing_km_vals, tau->arr_size);
    set_var(&py_tau->tau_threshold_vals, &tau->tau_threshold_vals, tau->arr_size);

    if (tau->T_fwd == NULL)
    {
        tmp = py_tau->p_norm_fwd_vals;
        Py_INCREF(Py_None);
        py_tau->p_norm_fwd_vals = Py_None;
        Py_XDECREF(tmp);

        tmp = py_tau->phase_fwd_vals;
        Py_INCREF(Py_None);
        py_tau->phase_fwd_vals = Py_None;
        Py_XDECREF(tmp);

        tmp = py_tau->tau_fwd_vals;
        Py_INCREF(Py_None);
        py_tau->tau_fwd_vals = Py_None;
        Py_XDECREF(tmp);
    }
    else
    {
        set_var(&py_tau->p_norm_fwd_vals, &tau->p_norm_fwd_vals, tau->arr_size);
        set_var(&py_tau->phase_fwd_vals, &tau->phase_fwd_vals, tau->arr_size);
        set_var(&py_tau->tau_fwd_vals, &tau->tau_fwd_vals, tau->arr_size);
    }
}

/*  To edit:
    PyObject          *rev_info;
    PyObject          *history;
*/

static void
rssringoccs_Get_Py_Perturb(rssringoccs_TAUObj *tau, PyObject *perturb)
{
    PyObject *iter;
    PyObject *next;
    unsigned int n;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    /*  Check that the input perturb is a list with 5 elements.               */
    if (!perturb)
    {
        for (n=0; n<5; ++n)
            tau->perturb[n] = 0.0;
    }

    /*  If the user supplied a perturb list, parse it and extract values.     */
    else if (PyList_Check(perturb))
    {
        /*  If the list is not the correct size, raise an error.              */
        if (PyList_Size(perturb) != 5)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Get_Py_Perturb\n\n"
                "\rInput perturb is a list but does not have 5 entries.\n"
                "\rperturb must be a list of five real numbers.\n"
            );
            return;
        }

        iter = PyObject_GetIter(perturb);

        /*  Loop over the elements of the list, see if they can be converted  *
         *  to doubles, and store them in the tau->perturb variable.          */
        for (n = 0; n < 5; ++n)
        {
            next = PyIter_Next(iter);

            /*  If the element is an integer, convert to double and save it.  */
            if (PyLong_Check(next))
                tau->perturb[n] = PyLong_AsDouble(next);

            /*  Convert from Python float to C double with PyFloat_AsDouble.  */
            else if (PyFloat_Check(next))
                tau->perturb[n] = PyFloat_AsDouble(next);

            /*  Invalid data type for one of the entries. Return with error.  */
            else
            {
                tau->error_occurred = tmpl_True;
                tau->error_message = tmpl_strdup(
                    "\rError Encountered: rss_ringoccs\n"
                    "\r\trssringoccs_Get_Py_Perturb\n\n"
                    "\rInput perturb has entries that are not real numbers.\n"
                    "\rAll five entries for the perturb list must be numbers.\n"
                );
                return;
            }
        }
    }

    /*  The input was not a list. Return with error.                          */
    else
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Get_Py_Perturb\n\n"
            "\rInput perturb is not a list.\n"
        );
        return;
    }
}

static void rssringoccs_Get_Py_Range(rssringoccs_TAUObj *tau, PyObject *rngreq)
{
    PyObject *iter;
    PyObject *next;
    unsigned int n;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (rngreq == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_Tau_to_Py_Tau\n\n"
            "\rInput rngreq is NULL. Aborting.n"
        );
        return;
    }

    /*  If the rng variable is a string, make sure it is a legal value and    *
     *  try to extract the corresponding values in kilometers.                */
    if PyBytes_Check(rngreq)
        rssringoccs_Tau_Set_Range_From_String(PyBytes_AsString(rngreq), tau);

    /*  If the rng variable is a unicode object (type of string from python)  *
     *  make sure it is a legal value and try to extract the corresponding    *
     *  values in kilometers.                                                 */
    else if PyUnicode_Check(rngreq)

        /*  Convert the Python string to a C string via PyUnicode_AsUTF8. The *
         *  C API recommends not altering the string, so we create a copy of  *
         *  it using strcpy (from string.h).                                  */
        rssringoccs_Tau_Set_Range_From_String(PyUnicode_AsUTF8(rngreq), tau);

    /*  If the requested range is a list, try to parse the elements.          */
    else if PyList_Check(rngreq)
    {
        if (PyList_Size(rngreq) != 2)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Get_Py_Range\n\n"
                "\rInput range is a list but does not have 2 entries.\n"
                "\rrng must be a list of two real numbers.\n"
            );
            return;
        }

        iter = PyObject_GetIter(rngreq);

        for (n = 0; n < 2; ++n)
        {
            next = PyIter_Next(iter);

            /*  Try to parse the elements. Return with error if this fails.   */
            if (PyLong_Check(next))
                tau->rng_list[n] = PyLong_AsDouble(next);
            else if (PyFloat_Check(next))
                tau->rng_list[n] = PyFloat_AsDouble(next);
            else
            {
                tau->error_occurred = tmpl_True;
                tau->error_message = tmpl_strdup(
                    "\rError Encountered: rss_ringoccs\n"
                    "\r\trssringoccs_Get_Py_Range\n\n"
                    "\rInput rng has entries that are not real numbers.\n"
                    "\rBoth entries for the rng list must be numbers.\n"
                );
                return;
            }
        }
    }

    /*  Illegal rng requested. Return with error.                             */
    else
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Get_Py_Range\n\n"
            "\rrng must be a list of two real numbers or a string.\n"
            "\rAllowed strings are:\n"
            "\r\tall               [1.0, 400000.0]\n"
            "\r\tbesselbarnard     [120210.0, 120330.0]\n"
            "\r\tbessel-barnard    [120210.0, 120330.0]\n"
            "\r\tcringripples      [77690.0, 77760.0]\n"
            "\r\tencke             [132900.0, 134200.0]\n"
            "\r\tenckegap          [132900.0, 134200.0]\n"
            "\r\therschel          [118100.0, 118380.0]\n"
            "\r\therschelgap       [118100.0, 118380.0]\n"
            "\r\thuygens           [117650.0, 117950.0]\n"
            "\r\thuygensringlet    [117650.0, 117950.0]\n"
            "\r\tjanusepimetheus   [96200.0, 96800.0]\n"
            "\r\tjeffreys          [118900.0, 119000.0]\n"
            "\r\tjeffreysgap       [118900.0, 119000.0]\n"
            "\r\tkuiper            [119300.0, 119500.0]\n"
            "\r\tkuipergap         [119300.0, 119500.0]\n"
            "\r\tmaxwell           [87410.0, 87610.0]\n"
            "\r\tmaxwellringlet    [87410.0, 87610.0]\n"
            "\r\trussell           [118550.0, 118660.0]\n"
            "\r\trussellgap        [118550.0, 118660.0]\n"
            "\r\ttitan             [77870.0, 77930.0]\n"
            "\r\ttitanringlet      [77870.0, 77930.0\n\n"
        );
        return;
    }
}

static void
rssringoccs_Get_Py_Vars_From_Self(rssringoccs_TAUObj *tau, PyDiffrecObj *self)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (self == NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_C_Tau_to_Py_Tau\n\n"
            "\rInput self is NULL. Aborting.n"
        );
        return;
    }

    tau->sigma    = self->sigma;
    tau->bfac     = self->bfac;
    tau->ecc      = self->ecc;
    tau->peri     = self->peri;
    tau->use_fwd  = self->use_fwd;
    tau->use_norm = self->use_norm;
    tau->verbose  = self->verbose;
}

/*  Deallocating function for the DiffractionCorrection class.                */
static void Diffrec_dealloc(PyDiffrecObj *self)
{
    Py_XDECREF(self->B_rad_vals);
    Py_XDECREF(self->D_km_vals);
    Py_XDECREF(self->F_km_vals);
    Py_XDECREF(self->f_sky_hz_vals);
    Py_XDECREF(self->p_norm_fwd_vals);
    Py_XDECREF(self->p_norm_vals);
    Py_XDECREF(self->phase_fwd_vals);
    Py_XDECREF(self->phase_rad_vals);
    Py_XDECREF(self->phase_vals);
    Py_XDECREF(self->phi_rad_vals);
    Py_XDECREF(self->phi_rl_rad_vals);
    Py_XDECREF(self->power_vals);
    Py_XDECREF(self->raw_tau_threshold_vals);
    Py_XDECREF(self->rev_info);
    Py_XDECREF(self->rho_corr_pole_km_vals);
    Py_XDECREF(self->rho_corr_timing_km_vals);
    Py_XDECREF(self->rho_dot_kms_vals);
    Py_XDECREF(self->rho_km_vals);
    Py_XDECREF(self->t_oet_spm_vals);
    Py_XDECREF(self->t_ret_spm_vals);
    Py_XDECREF(self->t_set_spm_vals);
    Py_XDECREF(self->tau_threshold_vals);
    Py_XDECREF(self->tau_vals);
    Py_XDECREF(self->tau_fwd_vals);
    Py_XDECREF(self->w_km_vals);
    Py_XDECREF(self->history);
    Py_XDECREF(self->input_vars);
    Py_XDECREF(self->input_kwds);
    Py_XDECREF(self->rx_km_vals);
    Py_XDECREF(self->ry_km_vals);
    Py_XDECREF(self->rz_km_vals);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

/*  The init function for the dirrection correction class. This is the        *
 *  equivalent of the __init__ function defined in a normal python class.     */
static int Diffrec_init(PyDiffrecObj *self, PyObject *args, PyObject *kwds)
{
    /*  Declare variables for a DLP and Tau object.                           */
    rssringoccs_DLPObj *dlp;
    rssringoccs_TAUObj *tau;

    /*  The list of the keywords accepted by the DiffractionCorrection class. *
     *  dlp and res are REQUIRED inputs, the rest are optional. If the user   *
     *  does not provide these optional keywords, we must set them ourselves. */
    static char *kwlist[] = {
        "dlp",
        "res",
        "rng",
        "wtype",
        "use_fwd",
        "use_norm",
        "verbose",
        "bfac",
        "sigma",
        "psitype",
        "res_factor",
        "ecc",
        "peri",
        "perturb",
        NULL
    };

    /*  Python objects needed throughout the computation.                     */
    PyObject *DLPInst;
    PyObject *tmp;
    PyObject *dlp_tmp;
    PyObject *rngreq;
    PyObject *perturb;

    /*  Set the default keyword options.                                      */

    /*  The kbmd20 is a new window, a modifed Kaiser-Bessel with alpha set to *
     *  two pi. The modification ensures the window goes to zero at its edges *
     *  while evaluating to one at the center, unlike the actual              *
     *  Kaiser-Bessel which is discontinuous at the edge of the window. The   *
     *  two pi factor is mostly guess work since it accurately reproduces the *
     *  PDS results. The real window used for that data is not known too me.  *
     *  The actual code for the window functions is in special_functions/     */
    self->wtype = "kbmd20";

    /*  Fresnel 4 is a new option, not mentioned in any of the papers but     *
     *  documented in our accompanying PDF. It uses Legendre polynomials to   *
     *  approximate the Fresnel kernel. It essentially takes Fresnels         *
     *  quadratic method to the next step, a quartic, hence the name. It is   *
     *  extremely fast (all of Rev007 takes less than a second) and very      *
     *  accurate for all but the most extreme occultations (like Rev133).     */
    self->psitype = "fresnel4";

    /*  Default range is "all", denoting [1.0, 400000.0]. We'll set later.    */
    rngreq = PyUnicode_FromString("all");

    /*  By default, forward computations are not run, FFTs are not used, and  *
     *  the run is silent (verbose is off).                                   */
    self->use_fwd = tmpl_False;
    self->verbose = tmpl_False;

    /*  Using the bfac guarantees accurate window sizes in the case of a poor *
     *  Allen deviation. Window normalization is also recommended since the   *
     *  integral is scaled by the width of the window, and hence for small    *
     *  window sizes the result might return close to zero.                   */
    self->bfac = tmpl_True;
    self->use_norm = tmpl_True;

    /*  The default sigma value is the one for Cassini.                       */
    self->sigma = 2.0e-13;

    /*  If res_factor was not set, set to 0.75. This value was specified by   *
     *  Essam Marouf as necessary to ensure the reconstruction matches the    *
     *  PDS results. No justification is known to me.                         */
    self->res_factor = 0.75;

    /*  The default geometry assumes the rings are circular, so we set both   *
     *  the eccentricity and the periapse to zero.                            */
    self->ecc = 0.0;
    self->peri = 0.0;

    /*  Default polynomial perturbation is off.                               */
    perturb = NULL;

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
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Od$OsppppdsdddO:", kwlist,
                                     &DLPInst,          &self->input_res,
                                     &rngreq,           &self->wtype,
                                     &self->use_fwd,    &self->use_norm,
                                     &self->verbose,    &self->bfac,
                                     &self->sigma,      &self->psitype,
                                     &self->res_factor, &self->ecc,
                                     &self->peri,       &perturb))
    {
        PyErr_Format(
            PyExc_TypeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rCould not parse input variables.\n\n"
            "\rInputs:\n"
            "\r\tDLPInst:  \tAn instance of the DLP Class.\n"
            "\r\tres:      \tRequested resolution in km (float).\n\n"
            "\rKeywords:\n"
            "\r\trng       \tThe requested range (str or list).\n"
            "\r\twtype     \tThe requested window type (str).\n"
            "\r\tuse_fwd   \tForward computation (bool).\n"
            "\r\tuse_norm  \tWindow normalization (bool).\n"
            "\r\tverbose   \tPrint status updates (bool).\n"
            "\r\tbfac      \tUse b-factor in window width (bool).\n"
            "\r\tsigma     \tThe Allen deviation (float).\n"
            "\r\tpsitype   \tRequested Frensel kernel approximation (str).\n"
            "\r\tres_factor\tScaling factor for resolution (float).\n"
            "\r\tecc       \tEccentricity of rings (bool).\n"
            "\r\tperi      \tPeriapse of rings (bool).\n"
            "\r\tperturb   \tRequested perturbation to Fresnel kernel (list).\n"
        );
        return -1;
    }

    if (self->verbose)
    {
        puts("Diffraction Correction:");
        puts("\tDiffraction Correction: Retrieving history from DLP...");
    }

    if (!PyObject_HasAttrString(DLPInst, "rev_info"))
    {
        PyErr_Format(
            PyExc_AttributeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rInput DLP Instance is missing the following attribute:\n"
            "\r\trev_info\n\n"
        );
        return -1;
    }
    else
        dlp_tmp = PyObject_GetAttrString(DLPInst, "rev_info");

    /*  Store the dlp history inside of the DiffractionCorrection class. This *
     *  tmp, Py_INCREF, Py_XDECREF method is recommended in the Python C-API  *
     *  documentation as a means of safely storing the variable.              */
    tmp = self->rev_info;
    Py_INCREF(dlp_tmp);
    self->rev_info = dlp_tmp;
    Py_XDECREF(tmp);

    /*  If verbose was set, print a status update.                            */
    if (self->verbose)
        puts("\tDiffraction Correction: Converting Py DLP to C DLP...");

    dlp = rssringoccs_Py_DLP_to_C_DLP(DLPInst);

    if (dlp == NULL)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rFailed to pass variables to C. rssringoccs_Py_DLP_to_C_DLP\n"
            "\rreturned NULL. Returning.\n\n"
        );
        return -1;
    }

    if (dlp->error_occurred)
    {
        if (dlp->error_message == NULL)
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
        {
            PyErr_Format(PyExc_RuntimeError, "%s", dlp->error_message);
            free(dlp->error_message);
        }
        free(dlp);
        return -1;
    }

    /*  If verbose was set, print a status update.                            */
    if (self->verbose)
        puts("\tDiffraction Correction: Creating C Tau object...");

    tau = rssringoccs_Create_TAUObj(dlp, self->input_res * self->res_factor);

    if (self->verbose)
        puts("\tDiffraction Correction: Passing Py variables to tau...");

    rssringoccs_Get_Py_Vars_From_Self(tau, self);
    rssringoccs_Get_Py_Perturb(tau, perturb);
    rssringoccs_Get_Py_Range(tau, rngreq);
    rssringoccs_Tau_Set_WType(self->wtype, tau);
    rssringoccs_Tau_Set_Psitype(self->psitype, tau);

    if (self->verbose)
        puts("\tDiffraction Correction: Running reconstruction...");

    rssringoccs_Reconstruction(tau);

    if (self->verbose)
        puts("\tDiffraction Correction: Converting C tau to Py tau...");

    rssringoccs_C_Tau_to_Py_Tau(self, tau);

    if (tau == NULL)
    {
        PyErr_Format(
            PyExc_RuntimeError,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tdiffrec.DiffractionCorrection\n\n"
            "\rrssringoccs_Create_TAUObj returned NULL for tau. Returning.\n\n"
        );

        return -1;
    }

    if (tau->error_occurred)
    {
        if (tau->error_message == NULL)
            PyErr_Format(
                PyExc_RuntimeError,
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\tdiffrec.DiffractionCorrection\n\n"
                "\rtau returned with error_occurred set to true but no\n"
                "\rerror message. Returning.\n\n"
            );
        else
            PyErr_Format(PyExc_RuntimeError, "%s\n", tau->error_message);

        rssringoccs_Destroy_Tau(&tau);
        return -1;
    }

    free(tau->psitype);
    free(tau->wtype);

    /*  We are now freeing the C tau object. The data pointers are still      *
     *  accessible via the self PyObject. Note, we are freeing the pointer to *
     *  the rssringoccs_TAUObj and NOT the pointers inside the object. The    *
     *  data is still available in self.                                      */
    free(tau);

    /*  Similarly, we free the DLP. This does not free the data from the      *
     *  input DLP PyObject. Those are also still available.                   */
    free(dlp);

    if (self->verbose)
        puts("\tDiffraction Correction: Building arguments dictionary...");

    dlp_tmp = Py_BuildValue(
        "{s:O,s:d}",
        "dlp_inst", PyObject_GetAttrString(DLPInst, "history"),
        "res",      self->input_res
    );

    tmp = self->input_vars;
    Py_INCREF(dlp_tmp);
    self->input_vars = dlp_tmp;
    Py_XDECREF(tmp);

    if (self->verbose)
        puts("\tDiffraction Correction: Building keywords dictionary...");

    dlp_tmp = Py_BuildValue(
        "{s:O,s:s,s:s,s:d,s:d,s:d,s:d,s:O,s:O}",
        "rng",        rngreq,
        "wtype",      self->wtype,
        "psitype",    self->psitype,
        "sigma",      self->sigma,
        "ecc",        self->ecc,
        "peri",       self->peri,
        "res_factor", self->res_factor,
        "use_norm",   PyBool_FromLong(self->use_norm),
        "bfac",       PyBool_FromLong(self->bfac)
    );

    tmp = self->input_kwds;
    Py_INCREF(dlp_tmp);
    self->input_kwds = dlp_tmp;
    Py_XDECREF(tmp);

    return 1;
}

static PyMemberDef Custom_members[] = {
    {
        "rho_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, rho_km_vals), 0,
        "Ring radius."
    },
    {
        "phase_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phase_rad_vals),
        0, "Raw diffracted phase."
    },
    {
        "B_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, B_rad_vals), 0,
        "Ring inclination angle."
    },
    {
        "D_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, D_km_vals), 0,
        "Spacecraft to ring-intercept point distance."
    },
    {
        "rx_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, rx_km_vals), 0,
        "x coordinate of the spacecraft in planetocentric frame."
    },
    {
        "ry_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, ry_km_vals), 0,
        "y coordinate of the spacecraft in planetocentric frame."
    },
    {
        "rz_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, rz_km_vals), 0,
        "z coordinate of the spacecraft in planetocentric frame."
    },
    {
        "F_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, F_km_vals), 0,
        "Fresnel scale."
    },
    {
        "f_sky_hz_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, f_sky_hz_vals), 0,
        "Frequency of the input signal"
    },
    {
        "p_norm_fwd_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, p_norm_fwd_vals),
        0, "Forward modeling of power"
    },
    {
        "p_norm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, p_norm_vals), 0,
        "Raw power data"
    },
    {
        "phase_fwd_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phase_fwd_vals),
        0, "Forward modeling of phase"
    },
    {
        "phase_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phase_vals), 0,
        "Diffraction corrected phase"
    },
    {
        "phi_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phi_rad_vals), 0,
        "Ring azimuth angle"
    },
    {
        "phi_rl_rad_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, phi_rl_rad_vals),
        0, "Ring longitude angle"
    },
    {
        "power_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, power_vals), 0,
        "Diffraction corrected power"
    },
    {
        "raw_tau_threshold_vals", T_OBJECT_EX,
        offsetof(PyDiffrecObj, raw_tau_threshold_vals), 0,
        "Raw threshold optical depth"
    },
    {
        "rev_info", T_OBJECT_EX, offsetof(PyDiffrecObj, rev_info), 0,
        "Information about the occultation"
    },
    {
        "rho_corr_pole_km_vals", T_OBJECT_EX,
        offsetof(PyDiffrecObj, rho_corr_pole_km_vals), 0,
        "Ring radius with pole correction."
    },
    {
        "rho_corr_timing_km_vals", T_OBJECT_EX,
        offsetof(PyDiffrecObj, rho_corr_timing_km_vals), 0,
        "Ring radius with timing correction."
    },
    {
        "rho_dot_kms_vals", T_OBJECT_EX,
        offsetof(PyDiffrecObj, rho_dot_kms_vals), 0,
        "Time derivative of the ring radius."
    },
    {
        "t_oet_spm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, t_oet_spm_vals),
        0, "Observed event time in seconds past midnight"
    },
    {
        "t_ret_spm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, t_ret_spm_vals),
        0, "Ring event time in seconds past midnight"
    },
    {
        "t_set_spm_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, t_set_spm_vals),
        0, "Spacecraft event time in seconds past midnight"
    },
    {
        "tau_threshold_vals", T_OBJECT_EX,
        offsetof(PyDiffrecObj, tau_threshold_vals), 0,
        "Diffraction corrected threshold optical depth"
    },
    {
        "tau_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, tau_vals), 0,
        "Optical depth"
    },
    {
        "tau_fwd_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, tau_fwd_vals), 0,
        "Optical depth"
    },
    {
        "w_km_vals", T_OBJECT_EX, offsetof(PyDiffrecObj, w_km_vals), 0,
        "window width as a function of ring radius"
    },
    {
        "history", T_OBJECT_EX, offsetof(PyDiffrecObj, history), 0,
        "History of the tau instance"
    },
    {
        "input_vars", T_OBJECT_EX, offsetof(PyDiffrecObj, input_vars), 0,
        "Dictionary of input arguments used to create this instance."
    },
    {
        "input_kwds", T_OBJECT_EX, offsetof(PyDiffrecObj, input_kwds), 0,
        "Dictionary of input keywords used to create this instance."
    },
    {
        "bfac", T_BOOL, offsetof(PyDiffrecObj, bfac), 0,
        "Use of b-factor in window width"
    },
    {
        "verbose", T_BOOL, offsetof(PyDiffrecObj, verbose), 0,
        "Print status updates"
    },
    {
        "use_norm", T_BOOL, offsetof(PyDiffrecObj, use_norm), 0,
        "Use of window normalization"
    },
    {
        "use_fwd", T_BOOL, offsetof(PyDiffrecObj, use_fwd), 0,
        "Forward modeling Boolean"
    },
    {
        "ecc", T_DOUBLE, offsetof(PyDiffrecObj, ecc), 0,
        "Eccentricity of Rings"
    },
    {
        "peri", T_DOUBLE, offsetof(PyDiffrecObj, peri), 0,
        "Periapse of Rings, azimuth angle in radians."
    },
    {
        "input_res",
        T_DOUBLE,
        offsetof(PyDiffrecObj, input_res),
        0,
        "User requested input resolution."
    },
    {
        "res_factor",
        T_DOUBLE,
        offsetof(PyDiffrecObj, res_factor),
        0,
        "User requested scale factor for the input resolution."
    },
    {
        "sigma",
        T_DOUBLE,
        offsetof(PyDiffrecObj, sigma),
        0,
        "Allen deviation."
    },
    {
        NULL
    }  /* Sentinel */
};

static PyObject *
DiffractionCorrection(PyDiffrecObj *self, PyObject *Py_UNUSED(ignored))
{
    return PyUnicode_FromFormat("DiffractionCorrection");
}

static PyMethodDef DiffractionCorrection_methods[] =
{
    {
        "DiffractionCorrection",
        (PyCFunction) DiffractionCorrection,
        METH_NOARGS,
        "Diffraction correction class."
    },
    {NULL}
};

static PyTypeObject DiffrecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "diffrec.DiffractionCorrection",
    .tp_doc = "Diffraction Correction class.",
    .tp_basicsize = sizeof(PyDiffrecObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) Diffrec_init,
    .tp_dealloc = (destructor) Diffrec_dealloc,
    .tp_members = Custom_members,
    .tp_methods = DiffractionCorrection_methods,
};

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "custom",
    .m_doc = "Module containing the rss_ringoccs class.",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_diffrec(void)
{
    PyObject *m;
    int pymod_addobj;
    if (PyType_Ready(&DiffrecType) < 0)
        return NULL;

    m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    Py_INCREF(&DiffrecType);
    pymod_addobj = PyModule_AddObject(m, "DiffractionCorrection",
                                      (PyObject *) &DiffrecType);

    if (pymod_addobj < 0)
    {
        Py_DECREF(&DiffrecType);
        Py_DECREF(m);
        return NULL;
    }

    import_array();
    return m;
}
