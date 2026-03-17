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
 *      Free all of the pointers in a Cal CSV object.                         *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CalCSV_Destroy_Members                                    *
 *  Purpose:                                                                  *
 *      Free's the memory in a Cal CSV and sets the members to NULL.          *
 *  Arguments:                                                                *
 *      cal (rssringoccs_CalCSV * const):                                     *
 *          A pointer to the Calibration CSV object we are destroying.        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      stdlib.h:                                                             *
 *          free:                                                             *
 *              free's memory allocated by malloc, calloc, or realloc.        *
 *  Method:                                                                   *
 *      Use the TMPL_FREE macro to free the data and set the pointers to NULL.*
 *  Notes:                                                                    *
 *      1.) This function checks for NULL pointers. If cal = NULL, nothing is *
 *          done. Members of the Cal CSV object that are NULL are skipped.    *
 *                                                                            *
 *      2.) To prevent double free's, the members of the Cal CSV object are   *
 *          set to NULL after free'ing.                                       *
 *                                                                            *
 *      3.) error_message is declared const char *, it is never free'd.       *
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

/*  Provides the TMPL_FREE macro for freeing a pointer and setting it to NULL.*/
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_CalCSV typedef found here.                                    */
#include <rss_ringoccs/include/types/rss_ringoccs_calcsv.h>

/*  Forward declaration / function prototype.                                 */
extern void rssringoccs_CalCSV_Destroy_Members(rssringoccs_CalCSV * const cal);

/*  Free's all members of an rssringoccs_CalCSV object except the             *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_CalCSV_Destroy_Members(rssringoccs_CalCSV * const cal)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!cal)
        return;

    /*  Free every pointer except the error message (which is declared const).*/
    TMPL_FREE(cal->history);
    TMPL_FREE(cal->t_oet_spm_vals);
    TMPL_FREE(cal->f_sky_pred_vals);
    TMPL_FREE(cal->f_sky_resid_fit_vals);
    TMPL_FREE(cal->p_free_vals);
}
/*  End of rssringoccs_CalCSV_Destroy_Members.                                */
