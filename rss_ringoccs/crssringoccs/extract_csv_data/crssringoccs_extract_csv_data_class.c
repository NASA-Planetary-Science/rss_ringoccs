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
#include "../crssringoccs.h"

static PyMemberDef extractcsvdata_members[] = {
    {
        "rho_km_vals", T_OBJECT_EX, offsetof(PyCSVObj, rho_km_vals), 0,
        "Ring radius."
    },
    {
        "phase_deg_vals", T_OBJECT_EX, offsetof(PyCSVObj, phase_deg_vals),
        0, "Raw diffracted phase."
    },
    {
        "B_deg_vals", T_OBJECT_EX, offsetof(PyCSVObj, B_deg_vals), 0,
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
        "raw_tau_vals", T_OBJECT_EX, offsetof(PyCSVObj, raw_tau_vals), 0,
        "Raw optical depth"
    },
    {
        "phi_deg_vals", T_OBJECT_EX, offsetof(PyCSVObj, phi_deg_vals), 0,
        "Ring azimuth angle"
    },
    {
        "phi_rl_deg_vals", T_OBJECT_EX, offsetof(PyCSVObj, phi_rl_deg_vals),
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
        "tau_power", T_OBJECT_EX, offsetof(PyCSVObj, tau_power), 0,
        "Optical power"
    },
    {
        "tau_phase", T_OBJECT_EX, offsetof(PyCSVObj, tau_phase), 0,
        "Optical phase"
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

PyTypeObject ExtractCSVDataType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "csvtools.ExtractCSVData",
    .tp_doc = "ExtractCSVData class.",
    .tp_basicsize = sizeof(PyCSVObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) ExtractCSVData_init,
    .tp_dealloc = (destructor) ExtractCSVData_dealloc,
    .tp_members = extractcsvdata_members,
    .tp_methods = ExtractCSVData_methods,
};
