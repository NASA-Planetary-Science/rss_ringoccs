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
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       March 10, 2026                                                *
 ******************************************************************************/
#include "../crssringoccs.h"

PyGetSetDef crssringoccs_DiffractionCorrection_GetSetters[] =
{
    {
        "p_norm_vals",
        crssringoccs_DiffractionCorrection_Get_P_Norm_Vals,
        crssringoccs_DiffractionCorrection_Set_P_Norm_Vals,
        "Normalized diffracted power from the input DLP object.",
        NULL
    },
    {
        "power_vals",
        crssringoccs_DiffractionCorrection_Get_Power_Vals,
        crssringoccs_DiffractionCorrection_Set_Power_Vals,
        "Reconstructed power from the input data.",
        NULL
    },
    {
        "p_fwd_vals",
        crssringoccs_DiffractionCorrection_Get_P_Fwd_Vals,
        crssringoccs_DiffractionCorrection_Set_P_Fwd_Vals,
        "Forward model power computed from the reconstructed data.",
        NULL
    },
    {
        "phase_norm_deg_vals",
        crssringoccs_DiffractionCorrection_Get_Phase_Norm_Deg_Vals,
        crssringoccs_DiffractionCorrection_Set_Phase_Norm_Deg_Vals,
        "Normalized diffracted phase from the input DLP object.",
        NULL
    },
    {
        "phase_deg_vals",
        crssringoccs_DiffractionCorrection_Get_Phase_Deg_Vals,
        crssringoccs_DiffractionCorrection_Set_Phase_Deg_Vals,
        "Reconstructed phase computed from the reconstructed data.",
        NULL
    },
    {
        "phase_fwd_deg_vals",
        crssringoccs_DiffractionCorrection_Get_Phase_Fwd_Deg_Vals,
        crssringoccs_DiffractionCorrection_Set_Phase_Fwd_Deg_Vals,
        "Forward model phase computed from the reconstructed data.",
        NULL
    },
    {
        "tau_norm_vals",
        crssringoccs_DiffractionCorrection_Get_Tau_Norm_Vals,
        crssringoccs_DiffractionCorrection_Set_Tau_Norm_Vals,
        "Normalized diffracted optical depth from the input DLP object.",
        NULL
    },
    {
        "tau_vals",
        crssringoccs_DiffractionCorrection_Get_Tau_Vals,
        crssringoccs_DiffractionCorrection_Set_Tau_Vals,
        "Reconstructed optical depth computed from the reconstructed data.",
        NULL
    },
    {
        "tau_fwd_vals",
        crssringoccs_DiffractionCorrection_Get_Tau_Fwd_Vals,
        crssringoccs_DiffractionCorrection_Set_Tau_Fwd_Vals,
        "Forward model optical depth computed from the reconstructed data.",
        NULL
    },
    {
        NULL
    }
};
