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
#include "../crssringoccs.h"

static PyMemberDef diffrec_members[] = {
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

PyTypeObject DiffrecType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "diffrec.DiffractionCorrection",
    .tp_doc = "Diffraction Correction class.",
    .tp_basicsize = sizeof(PyDiffrecObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) Diffrec_init,
    .tp_dealloc = (destructor) Diffrec_dealloc,
    .tp_members = diffrec_members,
    .tp_methods = DiffractionCorrection_methods,
};
