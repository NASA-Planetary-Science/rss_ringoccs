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
 *                          rss_ringoccs_cal_csv_init                         *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Initialize all members of the Cal CSV object to their zero values.    *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CalCSV_Init                                               *
 *  Purpose:                                                                  *
 *      Initializes all members of a CalCSV object to their zero values.      *
 *  Arguments:                                                                *
 *      cal (rssringoccs_CalCSV * const):                                     *
 *          A pointer to the Calibration CSV object we are initializing.      *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      None.                                                                 *
 *  Method:                                                                   *
 *      Set each pointer to NULL, each Boolean to False, and each number to 0.*
 *  Notes:                                                                    *
 *      1.) If cal is NULL, nothing is done.                                  *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_free.h:                                                          *
 *          Header file providing the TMPL_FREE macro.                        *
 *  2.) rss_ringoccs_calcsv.h:                                                *
 *          Header file containing the rssringoccs_CalCSV typedef.            *
 *  3.) stddef.h:                                                             *
 *          Standard header file providing the NULL macro.                    *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2026/03/17: Ryan Maguire                                                  *
 *      Added docstring, cleaned up a bit.                                    *
 ******************************************************************************/

/*  Booleans (True and False) provided here.                                  */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_CalCSV typedef found here.                                    */
#include <rss_ringoccs/include/types/rss_ringoccs_calcsv.h>

/*  NULL is found here.                                                       */
#include <stddef.h>

/*  Forward declaration / function prototype.                                 */
extern void rssringoccs_CalCSV_Init(rssringoccs_CalCSV * const cal);

/* Sets all variables in a CAL CSV to their default values.                   */
void rssringoccs_CalCSV_Init(rssringoccs_CalCSV * const cal)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!cal)
        return;

    /*  No history yet. Set this to NULL.                                     */
    cal->history = NULL;

    /*  All arrays should be empty, for now.                                  */
    cal->t_oet_spm_vals = NULL;
    cal->f_sky_pred_vals = NULL;
    cal->f_sky_resid_fit_vals = NULL;
    cal->p_free_vals = NULL;
    cal->n_elements = 0;

    /*  And no error has occurred, yet.                                       */
    cal->error_message = NULL;
    cal->error_occurred = tmpl_False;
}
/*  End of rssringoccs_CalCSV_Init.                                           */
