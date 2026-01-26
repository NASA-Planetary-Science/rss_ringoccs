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

/*  Booleans, interpolation, math routines, and more.                         */
#include <libtmpl/include/tmpl.h>

/*  Typedefs for CSV structs and function prototype given here.               */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  malloc and free are found here.                                           */
#include <stdlib.h>

/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  strtok function provided here.                                            */
#include <string.h>

rssringoccs_TauCSV *
rssringoccs_TauCSV_Extract(const char *filename, tmpl_Bool use_deprecated)
{
    /*  Pointer to the Tau struct.                                            */
    rssringoccs_TauCSV *tau;

    /*  File pointer for the CSV file.                                        */
    FILE *fp;

    /*  Allocate memory for the tau data.                                     */
    tau = malloc(sizeof(*tau));

    /*  Check if malloc failed.                                               */
    if (!tau)
        return NULL;

    /*  Initialize the pointers in the tau struct to NULL. The function       *
     *  rssringoccs_Destroy_TauCSV_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    rssringoccs_TauCSV_Init(tau);

    /*  Init sets use_deprecated to False by default. Set the correct value.  */
    tau->use_deprecated = use_deprecated;

    /*  Try to open the input file.                                           */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (!fp)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_Get_Tau\n\n"
            "fopen returned NULL. Failed to open file for reading.\n"
            "It is likely the filename is incorrect or does not exist.\n";

        return tau;
    }

    /*  Run a sanity check on the CSV file. It should have 12 or 13 columns.  */
    rssringoccs_TauCSV_Check_Column_Count(tau, fp);

    /*  Count the number of lines in the CSV file and malloc enough data for  *
     *  the the arrays (one for each column).                                 */
    rssringoccs_TauCSV_Malloc(tau, fp);

    /*  Read the data from the file pointer in to the CSV struct.             */
    rssringoccs_TauCSV_Read_Data(tau, fp);

    /*  Close the file.                                                       */
    fclose(fp);
    return tau;
}
/*  End of rssringoccs_TauCSV_Extract.                                        */
