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
 *      Free all of the pointers in a Cal CSV object.                         *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 31, 2020                                             *
 ******************************************************************************/

/*  Macro for freeing a pointer and setting it to NULL.                       */
#include <libtmpl/include/compat/tmpl_free.h>

/*  rssringoccs_CalCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Free's all members of an rssringoccs_CalCSV pointer except the            *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_CalCSV_Destroy_Members(rssringoccs_CalCSV *cal)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!cal)
        return;

    /*  Destroy every variable except the error_message.                      */
    TMPL_FREE(cal->t_oet_spm_vals);
    TMPL_FREE(cal->f_sky_pred_vals);
    TMPL_FREE(cal->f_sky_resid_fit_vals);
    TMPL_FREE(cal->p_free_vals);
    TMPL_FREE(cal->history);
}
/*  End of rssringoccs_CalCSV_Destroy_Members.                                */
