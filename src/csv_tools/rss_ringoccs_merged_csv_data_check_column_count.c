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
 *      Checks a DLPM CSV file and ensures it has the right number of columns.*
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 30, 2024                                            *
 ******************************************************************************/

/*  libtmpl provides Booleans, string duplicate, and CSV reading tools.       */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_utility.h>

/*  rssringoccs_DLPCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for checking the number of columns in a DLPM CSV file.           */
void
rssringoccs_MergedCSVData_Check_Column_Count(rssringoccs_MergedCSVData *dlpm,
                                             FILE *fp)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    size_t column_count;

    /*  Nothing to do if the DLPM object is NULL.                             */
    if (!dlpm)
        return;

    /*  Similarly, if an error already occurred, abort.                       */
    if (dlpm->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        dlpm->error_occurred = tmpl_True;
        dlpm->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_MergedCSVData_Check_Column_Count\n\n"
            "Input file is NULL. Aborting.\n"
        );

        return;
    }

    /*  libtmpl has tools for counting columns. Use this.                     */
    column_count = tmpl_CSV_Column_Count(fp);

    /*  There should be 19 columns in a merged DLP CSV file.                  */
    if (column_count != 19)
    {
        dlpm->error_occurred = tmpl_True;
        dlpm->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_MergedCSVData_Check_Column_Count\n\n"
            "DLPM CSV file does not have 19 columns. Aborting.\n"
        );
    }
}
/*  End of rssringoccs_MergeCSVData_Check_Column_Count.                       */
