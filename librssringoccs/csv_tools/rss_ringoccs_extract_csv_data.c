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

#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

rssringoccs_CSVData *
rssringoccs_Extract_CSV_Data(const char *geo,
                             const char *cal,
                             const char *dlp,
                             const char *tau,
                             rssringoccs_Bool use_deprecated)
{
    rssringoccs_GeoCSV  *geo_dat;
    rssringoccs_DLPCSV  *dlp_dat;
    rssringoccs_CalCSV  *cal_dat;
    rssringoccs_TauCSV  *tau_dat;
    rssringoccs_CSVData *csv_data;

    csv_data = malloc(sizeof(*csv_data));

    if (csv_data == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Extract_CSV_Data\n\n"
             "Malloc failed and returned NULL for csv_data. Returning.\n");
        return NULL;
    }

    geo_dat = rssringoccs_Get_Geo(geo, use_deprecated);
    if (geo_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Geo returned NULL for geo_dat. Aborting.\n"
        );
        return csv_data;
    }

    dlp_dat = rssringoccs_Get_DLP(dlp, use_deprecated);
    if (dlp_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_DLP returned NULL for dlp_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        return csv_data;
    }

    cal_dat = rssringoccs_Get_Cal(cal);
    if (cal_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Cal returned NULL for cal_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        return csv_data;
    }

    tau_dat = rssringoccs_Get_Tau(tau, use_deprecated);
    if (tau_dat == NULL)
    {
        csv_data->error_occurred = rssringoccs_True;
        csv_data->error_message = rssringoccs_strdup(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Extract_CSV_Data\n\n"
            "rssringoccs_Get_Tau returned NULL for tau_dat. Aborting.\n"
        );
        rssringoccs_Destroy_GeoCSV(&geo_dat);
        rssringoccs_Destroy_DLPCSV(&dlp_dat);
        rssringoccs_Destroy_CalCSV(&cal_dat);
        return csv_data;
    }

    rssringoccs_Destroy_GeoCSV(&geo_dat);
    rssringoccs_Destroy_DLPCSV(&dlp_dat);
    rssringoccs_Destroy_CalCSV(&cal_dat);
    rssringoccs_Destroy_TauCSV(&tau_dat);

    return csv_data;

}
