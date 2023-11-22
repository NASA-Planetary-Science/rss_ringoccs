/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       January 5, 2020                                               *
 ******************************************************************************/

#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

void rssringoccs_Tau_Destroy(rssringoccs_TAUObj **tau)
{
    rssringoccs_TAUObj *tau_inst = *tau;

    if (tau_inst == NULL)
        return;

    rssringoccs_Tau_Destroy_Members(tau_inst);

    if (tau_inst->error_message != NULL)
    {
        free(tau_inst->error_message);
        tau_inst->error_message = NULL;
    }

    free(tau_inst);
    *tau = NULL;
    return;
}
/*  End of rssringoccs_Tau_Destroy.                                           */
