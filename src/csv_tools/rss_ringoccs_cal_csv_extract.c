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
 *                        rss_ringoccs_cal_csv_extract                        *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Extracts all data from a Calibration CSV.                             *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CalCSV_Extract                                            *
 *  Purpose:                                                                  *
 *      Initializes a CalCSV and reads data from a given CAL.TAB file.        *
 *  Arguments:                                                                *
 *      filename (const char * const):                                        *
 *          The path for the CAL.TAB file being read.                         *
 *  Output:                                                                   *
 *      cal (rssringoccs_CalCSV *):                                           *
 *          The Calibration object containing the data from the CAL.TAB file. *
 *  Called Functions:                                                         *
 *      src/csv_tools/                                                        *
 *          rssringoccs_CalCSV_Init:                                          *
 *              Sets all members of a CalCSV to their zero values.            *
 *          rssringoccs_CalCSV_Check_Column_Count:                            *
 *              Checks a CSV file ensuring it has exactly four columns.       *
 *          rssringoccs_CalCSV_Malloc:                                        *
 *              Allocates memory for each of the CalCSV members.              *
 *          rssringoccs_CalCSV_Read_Data:                                     *
 *              Reads data from a FILE and writes it into a CalCSV object.    *
 *      stdlib.h:                                                             *
 *          malloc:                                                           *
 *              Dynamically allocates memory.                                 *
 *      stdio.h:                                                              *
 *          fopen:                                                            *
 *              Opens a file (with read permissions).                         *
 *          fclose:                                                           *
 *              Closes a file.                                                *
 *  Method:                                                                   *
 *      Allocate memory for a CalCSV, initialize its members, and then        *
 *      allocate memory for the arrays in the object. Once done, read the     *
 *      data from a CAL.TAB file and write it into the arrays of the struct.  *
 *  Notes:                                                                    *
 *      1.) If malloc fails to allocate memory for the CalCSV, NULL is        *
 *          returned. Check for this after using this function.               *
 *                                                                            *
 *      2.) If malloc succeeds, then a CalCSV pointer is returned. If any     *
 *          other error occurs (fopen fails, malloc cannot allocate memory    *
 *          for the arrays, etc.), then the error_occurred Boolean is set to  *
 *          True and an error message is stored in the error_message member.  *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  2.) tmpl_malloc.h:                                                        *
 *          Header file providing the TMPL_MALLOC macro.                      *
 *  3.) rss_ringoccs_csv_tools.h:                                             *
 *          Header file containing the rssringoccs_CalCSV typedef.            *
 *  4.) stdio.h:                                                              *
 *          Standard header file providing the FILE type, fopen, and fclose.  *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2026/03/17: Ryan Maguire                                                  *
 *      Added docstring, cleaned up a bit.                                    *
 ******************************************************************************/

/*  Booleans (True and False) provided here.                                  */
#include <libtmpl/include/tmpl_bool.h>

/*  TMPL_MALLOC macro found here.                                             */
#include <libtmpl/include/compat/tmpl_malloc.h>

/*  Typedefs for CSV structs and function prototype given here.               */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  Function for extracting the data from a CAL.TAB file.                     */
rssringoccs_CalCSV *rssringoccs_CalCSV_Extract(const char * const filename)
{
    /*  Pointer to a CalCSV struct.                                           */
    rssringoccs_CalCSV *cal;

    /*  File object for the file we're reading.                               */
    FILE *fp;

    /*  Allocate memory for the CalCSV object.                                */
    cal = TMPL_MALLOC(rssringoccs_CalCSV, 1);

    /*  Check if malloc failed.                                               */
    if (!cal)
        return NULL;

    /*  Initialize the pointers in the cal struct to NULL. The function       *
     *  rssringoccs_CalCSV_Destroy_Members will check which members are NULL  *
     *  and attempt to free those that aren't. Freeing a pointer that wasn't  *
     *  malloc'd will crash the program, hence this initialization.           */
    rssringoccs_CalCSV_Init(cal);

    /*  Try to open the input file with read permissions.                     */
    fp = fopen(filename, "r");

    /*  If fopen returned NULL, the file likely does not exist. Return error. */
    if (!fp)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message =
            "\r\nError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_CalCSV_Extract\n\n"
            "\rfopen returned NULL. Failed to open file for reading.\n"
            "\rIt is likely the filename is incorrect or does not exist.\n\n";

        return cal;
    }

    /*  Run a sanity check on the CSV file. It should have 4 columns.         */
    rssringoccs_CalCSV_Check_Column_Count(cal, fp);

    /*  Count the number of lines in the CSV file and malloc enough data for  *
     *  the four arrays (one for each column).                                */
    rssringoccs_CalCSV_Malloc(cal, fp);

    /*  Read the data from the file pointer in to the CSV struct.             */
    rssringoccs_CalCSV_Read_Data(cal, fp);

    /*  We're done with the file, close it.                                   */
    fclose(fp);

    /*  Return the CalCSV pointer. The previous functions inspect and set the *
     *  error_occurred Boolean to true should anything go wrong. The user     *
     *  should check this before using the data.                              */
    return cal;
}
/*  End of rssringoccs_CalCSV_Extract.                                        */
