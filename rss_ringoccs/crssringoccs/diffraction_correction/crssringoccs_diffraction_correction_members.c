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

PyMemberDef crssringoccs_DiffractionCorrection_Members[] = {
    {
        "rho_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, rho_km_vals),
        0,
        "Ring radius."
    },
    {
        "B_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, B_deg_vals),
        0,
        "Ring inclination angle."
    },
    {
        "D_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, D_km_vals),
        0,
        "Spacecraft to ring-intercept point distance."
    },
    {
        "rx_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, rx_km_vals),
        0,
        "x coordinate of the spacecraft in planetocentric frame."
    },
    {
        "ry_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, ry_km_vals),
        0,
        "y coordinate of the spacecraft in planetocentric frame."
    },
    {
        "rz_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, rz_km_vals),
        0,
        "z coordinate of the spacecraft in planetocentric frame."
    },
    {
        "F_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, F_km_vals),
        0,
        "Fresnel scale."
    },
    {
        "phi_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, phi_deg_vals),
        0,
        "Ring azimuth angle"
    },
    {
        "phi_rl_deg_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, phi_rl_deg_vals),
        0,
        "Ring longitude angle"
    },
    {
        "T_in",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, T_in), 0,
        "Diffraction profile."
    },
    {
        "T_out",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, T_out),
        0,
        "Diffraction corrected profile."
    },
    {
        "T_fwd",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, T_fwd),
        0,
        "Forward modeling profile."
    },
    {
        "rho_corr_pole_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, rho_corr_pole_km_vals),
        0,
        "Ring radius with pole correction."
    },
    {
        "rho_corr_timing_km_vals", T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, rho_corr_timing_km_vals),
        0,
        "Ring radius with timing correction."
    },
    {
        "rho_dot_kms_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, rho_dot_kms_vals),
        0,
        "Time derivative of the ring radius."
    },
    {
        "t_oet_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, t_oet_spm_vals),
        0,
        "Observed event time in seconds past midnight"
    },
    {
        "t_ret_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, t_ret_spm_vals),
        0,
        "Ring event time in seconds past midnight"
    },
    {
        "t_set_spm_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, t_set_spm_vals),
        0,
        "Spacecraft event time in seconds past midnight"
    },
    {
        "tau_threshold_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, tau_threshold_vals),
        0,
        "Diffraction corrected threshold optical depth"
    },
    {
        "w_km_vals",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, w_km_vals),
        0,
        "window width as a function of ring radius"
    },
    {
        "outfiles",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, outfiles),
        0,
        "TAB files for the Tau object."
    },
    {
        "input_vars",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, input_vars),
        0,
        "Dictionary of input arguments used to create this instance."
    },
    {
        "input_kwds",
        T_OBJECT_EX,
        offsetof(crssringoccs_PyDiffrecObj, input_kwds),
        0,
        "Dictionary of input keywords used to create this instance."
    },
    {
        "bfac",
        T_BOOL,
        offsetof(crssringoccs_PyDiffrecObj, bfac),
        0,
        "Use of b-factor in window width"
    },
    {
        "verbose",
        T_BOOL,
        offsetof(crssringoccs_PyDiffrecObj, verbose),
        0,
        "Print status updates"
    },
    {
        "use_norm",
        T_BOOL,
        offsetof(crssringoccs_PyDiffrecObj, use_norm),
        0,
        "Use of window normalization"
    },
    {
        "use_fwd",
        T_BOOL,
        offsetof(crssringoccs_PyDiffrecObj, use_fwd),
        0,
        "Forward modeling Boolean"
    },
    {
        "ecc",
        T_DOUBLE,
        offsetof(crssringoccs_PyDiffrecObj, ecc),
        0,
        "Eccentricity of Rings"
    },
    {
        "peri",
        T_DOUBLE,
        offsetof(crssringoccs_PyDiffrecObj, peri),
        0,
        "Periapse of Rings, azimuth angle in radians."
    },
    {
        "input_res",
        T_DOUBLE,
        offsetof(crssringoccs_PyDiffrecObj, input_res),
        0,
        "User requested input resolution."
    },
    {
        "res_factor",
        T_DOUBLE,
        offsetof(crssringoccs_PyDiffrecObj, res_factor),
        0,
        "User requested scale factor for the input resolution."
    },
    {
        "sigma",
        T_DOUBLE,
        offsetof(crssringoccs_PyDiffrecObj, sigma),
        0,
        "Allen deviation."
    },
    {
        NULL
    }
};
