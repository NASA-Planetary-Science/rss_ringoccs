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
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_Destroy_CalCSV(rssringoccs_CalCSV **cal)
{
    rssringoccs_CalCSV *cal_inst = *cal;

    if (cal_inst == NULL)
        return;

    rssringoccs_Destroy_CalCSV_Members(cal_inst);

    if (cal_inst->error_message != NULL)
    {
        free(cal_inst->error_message);
        cal_inst->error_message = NULL;
    }

    free(cal_inst);
    *cal = NULL;
    return;
}
