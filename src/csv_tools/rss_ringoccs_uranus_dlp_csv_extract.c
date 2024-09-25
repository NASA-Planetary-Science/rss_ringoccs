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
 *      Extracts all diffraction data from a DLP CSV.                         *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 24, 2024                                            *
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

/*  Function for extracting the data from a DLP.TAB file.                     */
rssringoccs_UranusDLPCSV *
rssringoccs_UranusDLPCSV_Extract(const char *filename, tmpl_Bool in_radians)
{
    /*  Pointer to the DLP struct.                                            */
    rssringoccs_UranusDLPCSV *dlp;

    /*  File pointer for the CSV file.                                        */
    FILE *fp;

    /*  Allocate memory for the dlp data.                                     */
    dlp = malloc(sizeof(*dlp));

    /*  Check if malloc failed.                                               */
    if (!dlp)
        return NULL;

    /*  Initialize the pointers in the dlp struct to NULL. The function       *
     *  rssringoccs_UranusDLPCSV_Destroy_Members will check which members are *
     *  NULL and attempt to free those that aren't. Freeing a pointer that    *
     *  wasn't malloc'd will crash the program, hence this initialization.    */
    rssringoccs_UranusDLPCSV_Init(dlp);

    dlp->in_radians = in_radians;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (!fp)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusDLPCSV_Extract\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n"
        );

        return dlp;
    }

    /*  Run a sanity check on the CSV file. It should have 17 columns.        */
    rssringoccs_UranusDLPCSV_Check_Column_Count(dlp, fp);

    /*  Count the number of lines in the CSV file and malloc enough data for  *
     *  the the arrays (one for each column).                                 */
    rssringoccs_UranusDLPCSV_Malloc(dlp, fp);

    /*  Read the data from the file pointer in to the CSV struct.             */
    rssringoccs_UranusDLPCSV_Read_Data(dlp, fp);

    /*  Close the file.                                                       */
    fclose(fp);
    return dlp;
}
/*  End of rssringoccs_UranusDLPCSV_Extract.                                  */