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
 *      Function for free'ing a DLPM CSV object and free'ing all of the       *
 *      pointers contained inside the struct.                                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 30, 2024                                            *
 ******************************************************************************/

/*  Booleans (True and False) provided here.                                  */
#include <libtmpl/include/tmpl_bool.h>

/*  Provides the TMPL_FREE macro for freeing a pointer and setting it to NULL.*/
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_MergedCSVData typedef here, and function prototype given.     */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Function for freeing the memory in a MergedCSVData object.                */
void rssringoccs_MergedCSVData_Destroy(rssringoccs_MergedCSVData ** const dlpm)
{
    /*  If the input pointer is NULL, do not attempt to free it.              */
    if (!dlpm)
        return;

    /*  If this pointer is NULL, there is nothing to free.                    */
    if (!(*dlpm))
        return;

    /*  Free all of the members inside the DLPM object.                       */
    rssringoccs_MergedCSVData_Destroy_Members(*dlpm);

    /*  The error_message member is a pointer to a constant string, it does   *
     *  not need to be freed. Set the pointer to NULL to avoid reading it.    */
    (*dlpm)->error_occurred = tmpl_False;
    (*dlpm)->error_message = NULL;

    /*  Free the DLPCSV object pointer and set it to NULL to prevent trying   *
     *  to free the pointer twice.                                            */
    TMPL_FREE(*dlpm);
}
/*  End of rssringoccs_MergedCSVData_Destroy.                                 */
