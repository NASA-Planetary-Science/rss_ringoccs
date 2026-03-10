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
        NULL
    }
};
