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
 *                   rss_ringoccs_cal_csv_check_column_count                  *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Checks a Cal CSV file and ensures it has the right number of columns. *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CalCSV_Check_Column_Count                                 *
 *  Purpose:                                                                  *
 *      Checks the number of columns in a Calibration CSV file.               *
 *  Arguments:                                                                *
 *      cal (rssringoccs_CalCSV * TMPL_RESTRICT const):                       *
 *          A pointer to the Calibration CSV object we are inspecting.        *
 *      fp (FILE * TMPL_RESTRICT const):                                      *
 *          The file pointer for the CSV file.                                *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      tmpl_utility.h:                                                       *
 *          tmpl_CSV_Column_Count:                                            *
 *              Counts the number of columns in a comma-separated-file.       *
 *  Method:                                                                   *
 *      Check for basic errors, and then use libtmpl to count the number of   *
 *      columns in the input CSV file. Raise an error if this count is not 4, *
 *      corresponding to time, predicted frequency, sky frequency, and the    *
 *      free space power values, respectively.                                *
 *  Notes:                                                                    *
 *      1.) This function checks for NULL pointers and inspects the error     *
 *          Boolean. Nothing is done if cal = NULL or error_occurred = True.  *
 *                                                                            *
 *      2.) libtmpl's column count function rewinds the FILE back to the      *
 *          start. The user does not need to manually do this themselves.     *
 *                                                                            *
 *      3.) The error_occurred Boolean is set to True if either the input     *
 *          file is NULL, or if the column count is not 4. Inspect this after *
 *          calling this function.                                            *
 *                                                                            *
 *      4.) Both parameters are declared with the TMPL_RESTRICT macro. When   *
 *          compilers supporting C99 (or higher) are used, this expands to    *
 *          "restrict." Because of this, cal and fp must pointer to different *
 *          objects. This should be the case regardless to properly use this. *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_config.h:                                                        *
 *          Header providing the TMPL_RESTRICT macro.                         *
 *  2.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans.                                   *
 *  3.) tmpl_utility.h:                                                       *
 *          CSV tools, including column count, given here.                    *
 *  4.) rss_ringoccs_calcsv.h:                                                *
 *          Header file containing the rssringoccs_CalCSV typedef.            *
 *  5.) stdio.h:                                                              *
 *          Standard library header providing the FILE and size_t types.      *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2026/03/17: Ryan Maguire                                                  *
 *      Added docstring, cleaned up a bit.                                    *
 *  2026/03/26: Ryan Maguire                                                  *
 *      Added the TMPL_RESTRICT macro to the parameters.                      *
 ******************************************************************************/

/*  TMPL_RESTRICT macro found here.                                           */
#include <libtmpl/include/tmpl_config.h>

/*  Booleans (True and False) provided here.                                  */
#include <libtmpl/include/tmpl_bool.h>

/*  CSV reading tools given here.                                             */
#include <libtmpl/include/tmpl_utility.h>

/*  rssringoccs_CalCSV typedef found here.                                    */
#include <rss_ringoccs/include/types/rss_ringoccs_calcsv.h>

/*  FILE and size_t typedefs provided here.                                   */
#include <stdio.h>

/*  Forward declaration / function prototype.                                 */
extern void
rssringoccs_CalCSV_Check_Column_Count(
    rssringoccs_CalCSV * TMPL_RESTRICT const cal,
    FILE * TMPL_RESTRICT const fp
);

/*  Function for checking the number of columns in a Cal CSV file.            */
void
rssringoccs_CalCSV_Check_Column_Count(rssringoccs_CalCSV * const cal,
                                      FILE * const fp)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    size_t column_count;

    /*  Nothing to do if the Cal object is NULL.                              */
    if (!cal)
        return;

    /*  Similarly, if an error already occurred, abort.                       */
    if (cal->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CalCSV_Check_Column_Count\n\n"
            "Input file is NULL.\n\n";

        return;
    }

    /*  libtmpl has tools for counting columns. Use this.                     */
    column_count = tmpl_CSV_Column_Count(fp);

    /*  There should be 4 columns. Check this.                                */
    if (column_count != 4)
    {
        cal->error_occurred = tmpl_True;
        cal->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CalCSV_Check_Column_Count\n\n"
            "Input CSV does not have 4 columns.\n\n";
    }
}
/*  End of rssringoccs_CalCSV_Check_Column_Count.                             */
