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
 *  Purpose:                                                                  *
 *      Extracts all geometry data from a Geo CSV.                            *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  malloc found here.                                                        */
#include <stdlib.h>

/*  libtmpl provides Booleans.                                                */
#include <libtmpl/include/tmpl_bool.h>

/*  Typedefs for CSV structs and function prototype given here.               */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for extracting the data from a GEO.TAB file.                     */
rssringoccs_GeoCSV *
rssringoccs_GeoCSV_Extract(const char *filename, tmpl_Bool use_deprecated)
{
    /*  Pointer to the Geo struct.                                            */
    rssringoccs_GeoCSV *geo;

    /*  File pointer for the CSV file.                                        */
    FILE *fp;

    /*  Allocate memory for the geo data.                                     */
    geo = malloc(sizeof(*geo));

    /*  Check if malloc failed.                                               */
    if (!geo)
        return NULL;

    /*  Initialize the pointers in the geo struct to NULL. The function       *
     *  rssringoccs_GeoCSV_Destroy_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    rssringoccs_GeoCSV_Init(geo);

    /*  Init sets use_deprecated to False by default. Set the correct value.  */
    geo->use_deprecated = use_deprecated;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (!fp)
    {
        geo->error_occurred = tmpl_True;
        geo->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_GeoCSV_Extract\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n";

        return geo;
    }

    /*  Run a sanity check on the CSV file. It should have 18 or 19 columns.  */
    rssringoccs_GeoCSV_Check_Column_Count(geo, fp);

    /*  Count the number of lines in the CSV file and malloc enough data for  *
     *  the the arrays (one for each column).                                 */
    rssringoccs_GeoCSV_Malloc(geo, fp);

    /*  Read the data from the file pointer in to the CSV struct.             */
    rssringoccs_GeoCSV_Read_Data(geo, fp);

    /*  Close the file.                                                       */
    fclose(fp);
    return geo;
}
/*  End of rssringoccs_GeoCSV_Extract.                                        */
