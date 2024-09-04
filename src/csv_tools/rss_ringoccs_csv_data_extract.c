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
 *  Purpose:                                                                  *
 *      Extracts the data from CSV files to mimic a DLP object. This object   *
 *      can then be used for diffraction reconstruction.                      *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Booleans, interpolation, math routines, and more.                         */
#include <libtmpl/include/tmpl.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  free and malloc are found here.                                           */
#include <stdlib.h>

/*  Extracts all data from CSV (.TAB) files and creates a DLP-like struct.    */
rssringoccs_CSVData *
rssringoccs_CSVData_Extract(const char *geo,
                            const char *cal,
                            const char *dlp,
                            const char *tau,
                            tmpl_Bool use_deprecated)
{
    /*  A pointer to the CSV object.                                          */
    rssringoccs_CSVData *csv;

    /*  Allocate memory for the CSV object.                                   */
    csv = malloc(sizeof(*csv));

    /*  Check if malloc failed.                                               */
    if (!csv)
        return NULL;

    /*  Initialize the members to NULL. This will prevent functions from      *
     *  trying to free pointers that weren't malloc'd in the event of error.  */
    rssringoccs_CSVData_Init(csv);

    csv->use_deprecated = use_deprecated;

    rssringoccs_CSVData_Extract_Geo(csv, geo);
    rssringoccs_CSVData_Extract_Cal(csv, cal);
    rssringoccs_CSVData_Extract_DLP(csv, dlp);

    /*  Extract the data from the TAU.TAB file.                               */
    if (tau)
        rssringoccs_CSVData_Extract_Tau(csv, tau);

    rssringoccs_CSVData_Steal_DLP_Data(csv);
    rssringoccs_CSVData_Malloc(csv);
    rssringoccs_CSVData_Check_Chord_Occ(csv);
    rssringoccs_CSVData_Interpolate_Geo(csv);
    rssringoccs_CSVData_Interpolate_Cal(csv);

    /*  Interpolate the Tau data if requested.                                */
    if (tau)
        rssringoccs_CSVData_Interpolate_Tau(csv);

    /*  Free the Geo, Cal, and TAU structs. The CSV struct stole several      *
     *  pointers from the DLP struct, so do not destroy this.                 */
    rssringoccs_GeoCSV_Destroy(&(csv->geo));
    rssringoccs_CalCSV_Destroy(&(csv->cal));
    rssringoccs_TauCSV_Destroy(&(csv->tau));
    free(csv->dlp);
    csv->dlp = NULL;
    return csv;
}
