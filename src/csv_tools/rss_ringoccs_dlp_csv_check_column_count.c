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
 *      Checks a DLP CSV file and ensures it has the right number of columns. *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************/

/*  libtmpl provides Booleans and CSV reading tools.                          */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_utility.h>

/*  rssringoccs_DLPCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for checking the number of columns in a DLP CSV file.            */
void rssringoccs_DLPCSV_Check_Column_Count(rssringoccs_DLPCSV *dlp, FILE *fp)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    size_t column_count;

    /*  Nothing to do if the DLP object is NULL.                              */
    if (!dlp)
        return;

    /*  Similarly, if an error already occurred, abort.                       */
    if (dlp->error_occurred)
        return;

    /*  If the file pointer is NULL, the user likely called this function by  *
     *  mistake. Treat this as an error and abort.                            */
    if (!fp)
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_DLPCSV_Check_Column_Count\n\n"
            "Input file is NULL.\n";

        return;
    }

    /*  libtmpl has tools for counting columns. Use this.                     */
    column_count = tmpl_CSV_Column_Count(fp);

    /*  If use_deprecated was set to true, column_count must be 12. Check.    */
    if ((column_count != 12) && (dlp->use_deprecated))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_DLPCSV_Check_Column_Count\n\n"
            "use_deprecated is true but the CSV does not have 12 columns.\n";
    }

    /*  And if use_deprecated is false, we need 13 column. Check this.        */
    else if ((column_count != 13) && (!dlp->use_deprecated))
    {
        dlp->error_occurred = tmpl_True;
        dlp->error_message =
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_DLPCSV_Check_Column_Count\n\n"
            "use_deprecated is false but the CSV does not have 13 columns.\n";
    }
}
/*  End of rssringoccs_DLPCSV_Check_Column_Count.                             */
