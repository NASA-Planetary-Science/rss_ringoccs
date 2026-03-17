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
 *                    rss_ringoccs_cal_csv_destroy_members                    *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Function for free'ing a Cal CSV object and free'ing all of the        *
 *      pointers contained inside the struct.                                 *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CalCSV_Destroy                                            *
 *  Purpose:                                                                  *
 *      Destroys all of the data associated with a Calibration CSV.           *
 *  Arguments:                                                                *
 *      cal (rssringoccs_CalCSV ** const):                                    *
 *          A pointer to the Calibration CSV object we are destroying.        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      src/csv_tools/                                                        *
 *          rssringoccs_CalCSV_Destroy_Members:                               *
 *              Free's the individual members in a Calibration CSV object.    *
 *      stdlib.h:                                                             *
 *          free:                                                             *
 *              free's memory allocated by malloc, calloc, or realloc.        *
 *  Method:                                                                   *
 *      Use rssringoccs_CalCSV_Destroy_Members to free all of the members in  *
 *      the Cal CSV, and then free the Cal CSV object itself.                 *
 *  Notes:                                                                    *
 *      1.) This function checks for NULL pointers. If cal = NULL, nothing is *
 *          done. Members of the Cal CSV object that are NULL are skipped.    *
 *                                                                            *
 *      2.) To prevent double free's, the members of the Cal CSV object are   *
 *          set to NULL after free'ing. cal itself is set to NULL as well.    *
 *                                                                            *
 *      3.) error_message is declared const char *, it is never free'd.       *
 *          Instead, we simply set this pointer to NULL.                      *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_free.h:                                                          *
 *          Header file providing the TMPL_FREE macro.                        *
 *  3.) rss_ringoccs_calcsv.h:                                                *
 *          Header file containing the rssringoccs_CalCSV typedef.            *
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

/*  Provides the TMPL_FREE macro for freeing a pointer and setting it to NULL.*/
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_CalCSV typedef found here.                                    */
#include <rss_ringoccs/include/types/rss_ringoccs_calcsv.h>

/*  NULL macro defined here.                                                  */
#include <stddef.h>

/*  Function prototype / forward declaration.                                 */
extern void rssringoccs_CalCSV_Destroy(rssringoccs_CalCSV ** const cal);

/*  Tell the compiler about the main destructor function.                     */
extern void rssringoccs_CalCSV_Destroy_Members(rssringoccs_CalCSV * const cal);

/*  Function for freeing the memory in a CalCSV object.                       */
void rssringoccs_CalCSV_Destroy(rssringoccs_CalCSV ** const cal)
{
    /*  If the input pointer is NULL, simply return.                          */
    if (!cal)
        return;

    /*  If this is NULL, there's no need to free it. Return.                  */
    if (!(*cal))
        return;

    /*  Free all of the pointers inside the CalCSV object.                    */
    rssringoccs_CalCSV_Destroy_Members(*cal);

    /*  The error_message member is a pointer to a constant string, it does   *
     *  not need to be freed. Set the pointer to NULL to avoid reading it.    */
    (*cal)->error_occurred = tmpl_False;
    (*cal)->error_message = NULL;

    /*  Free the CalCSV pointer and set it to NULL to prevent freeing twice.  */
    TMPL_FREE(*cal);
}
/*  End of rssringoccs_CalCSV_Destroy.                                        */
