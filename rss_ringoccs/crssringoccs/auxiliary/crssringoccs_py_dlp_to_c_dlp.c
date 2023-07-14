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

/*  NULL is defined here.                                                     */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  tmpl_strdup function declared here.                                       */
#include <libtmpl/include/tmpl_string.h>

/*  Function prototype and typedefs for structs given here.                   */
#include "../crssringoccs.h"

/*  Avoid warnings about deprecated Numpy API versions.                       */
#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

/*  Numpy header files.                                                       */
#include <numpy/ndarraytypes.h>
#include <numpy/ufuncobject.h>

rssringoccs_DLPObj *crssringoccs_Py_DLP_to_C_DLP(PyObject *py_dlp)
{
    PyObject *tmp;
    PyObject *arr;
    rssringoccs_DLPObj *dlp;

    if (PyArray_API == NULL)
        import_array();

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
            "\r\tcrssringoccs_Py_DLP_to_C_DLP\n\n"
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
