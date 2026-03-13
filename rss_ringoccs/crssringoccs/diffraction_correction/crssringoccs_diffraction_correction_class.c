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

PyTypeObject crssringoccs_DiffractionCorrection = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "DiffractionCorrection",
    .tp_doc =
        "\r\tPurpose:\n"
        "\r\t\tPerforms diffraction correction on a diffraction limited\n"
        "\r\t\tprofile (DLP). This data can come from CSV files (using the\n"
        "\r\t\tCSV classes such as CassiniCSVData, VoyagerCSVData, or\n"
        "\r\t\tMergedCSVData), or from the DiffractionLimitedProfile class.\n"
        "\r\tArguments:\n"
        "\r\t\tdata (DiffractionLimitProfile):\n"
        "\r\t\t\tAny Python object that contains the following arrays:\n"
        "\r\t\t\t\tp_norm_vals:\n"
        "\r\t\t\t\t\tNormalized diffracted power.\n"
        "\r\t\t\t\tphase_deg_vals:\n"
        "\r\t\t\t\t\tPhase, in degrees, of the complex diffracted signal.\n"
        "\r\t\t\t\trho_km_vals:\n"
        "\r\t\t\t\t\tRadius of the ring intercept point, in kilometers.\n"
        "\r\t\t\t\tphi_deg_vals:\n"
        "\r\t\t\t\t\tAzimuth angle, in degrees, of the ring intercept point.\n"
        "\r\t\t\t\tphi_rl_deg_vals:\n"
        "\r\t\t\t\t\tRing longitude angle, in degrees.\n"
        "\r\t\t\t\tB_deg_vals:\n"
        "\r\t\t\t\t\tOpening angle, in degrees, of the ring plane.\n"
        "\r\t\t\t\tD_km_vals:\n"
        "\r\t\t\t\t\tSpacecraft to ring intercept point distance.\n"
        "\r\t\t\t\tf_sky_hz_vals:\n"
        "\r\t\t\t\t\tSky frequency of the diffracted radio wave, in Hertz.\n"
        "\r\t\t\t\trho_dot_kms_vals:\n"
        "\r\t\t\t\t\tRadial velocity of the ring intercept point (drho/dt).\n"
        "\r\t\t\t\tt_oet_spm_vals:\n"
        "\r\t\t\t\t\tObserved event time, in seconds past midnight.\n"
        "\r\t\t\t\tt_ret_spm_vals:\n"
        "\r\t\t\t\t\tRing event time, in seconds past midnight.\n"
        "\r\t\t\t\tt_set_spm_vals:\n"
        "\r\t\t\t\t\tSpacecraft event time, in seconds past midnight."
        "\r\t\t\t\trx_km_vals:\n"
        "\r\t\t\t\t\tx-coordinate of the spacecraft, in kilometers.\n"
        "\r\t\t\t\try_km_vals:\n"
        "\r\t\t\t\t\ty-coordinate of the spacecraft, in kilometers.\n"
        "\r\t\t\t\trz_km_vals:\n"
        "\r\t\t\t\t\tz-coordinate of the spacecraft, in kilometers.\n"
        "\r\t\t\t\trho_corr_pole_km_vals:\n"
        "\r\t\t\t\t\trho_km_vals with pole corrections.\n"
        "\r\t\t\t\trho_corr_timing_km_vals:\n"
        "\r\t\t\t\t\trho_km_vals with timing corrections.\n"
        "\r\t\t\t\traw_tau_threshold_vals:\n"
        "\r\t\t\t\t\tThreshold optical depth for the diffracted signal.\n"
        "\r\t\tresolution_km (float):"
        "\r\t\t\tThe requested resolution for the reconstruction.\n"
        "\r\tKeywords:\n",
    .tp_basicsize = sizeof(crssringoccs_PyDiffrecObj),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = crssringoccs_DiffractionCorrection_New,
    .tp_init = (initproc)crssringoccs_DiffractionCorrection_Init,
    .tp_dealloc = (destructor)crssringoccs_DiffractionCorrection_Destroy,
    .tp_getset = crssringoccs_DiffractionCorrection_GetSetters,
    .tp_members = crssringoccs_DiffractionCorrection_Members,
    .tp_methods = crssringoccs_DiffractionCorrection_Methods
};
