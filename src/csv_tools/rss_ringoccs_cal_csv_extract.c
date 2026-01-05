/******************************************************************************
 *                                 LICENSE                                    *
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
 *      Extracts all data from a Calibration CSV.                             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  malloc found here.                                                        */
#include <stdlib.h>

/*  libtmpl provides Booleans and string duplicate.                           */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>

/*  Typedefs for CSV structs and function prototype given here.               */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for extracting the data from a CAL.TAB file.                     */
rssringoccs_CalCSV *rssringoccs_CalCSV_Extract(const char *filename)
{
    /*  Pointer to a CalCSV struct.                                           */
    rssringoccs_CalCSV *cal;

    /*  File object for the file we're reading.                               */
    FILE *fp;

    /*  Allocate memory for the CalCSV object.                                */
    cal = malloc(sizeof(*cal));

    /*  Check if malloc failed.                                               */
    if (!cal)
        return NULL;

    /*  Initialize the pointers in the cal struct to NULL. The function       *
     *  rssringoccs_CalCSV_Destroy_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    rssringoccs_CalCSV_Init(cal);

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (!fp)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_CalCSV_Extract\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n";

        return cal;
    }

    /*  Run a sanity check on the CSV file. It should have 4 columns.         */
    rssringoccs_CalCSV_Check_Column_Count(cal, fp);

    /*  Count the number of lines in the CSV file and malloc enough data for  *
     *  the four arrays (one for each column).                                */
    rssringoccs_CalCSV_Malloc(cal, fp);

    /*  Read the data from the file pointer in to the CSV struct.             */
    rssringoccs_CalCSV_Read_Data(cal, fp);

    /*  Close the file and return the Cal object.                            */
    fclose(fp);
    return cal;
}
/*  End of rssringoccs_CalCSV_Extract.                                        */
