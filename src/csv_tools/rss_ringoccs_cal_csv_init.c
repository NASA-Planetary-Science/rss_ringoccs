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
 *      Initialize all members of the Cal CSV object to their zero values.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 1, 2024                                             *
 ******************************************************************************/

/*  NULL is found here.                                                       */
#include <stddef.h>

/*  Booleans provided here.                                                   */
#include <libtmpl/include/tmpl_bool.h>

/*  rssringoccs_CalCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/* Sets all variables in a CAL CSV to their default values.                   */
void rssringoccs_CalCSV_Init(rssringoccs_CalCSV *cal)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!cal)
        return;

    /*  All arrays should be empty, for now.                                  */
    cal->t_oet_spm_vals = NULL;
    cal->f_sky_pred_vals = NULL;
    cal->f_sky_resid_fit_vals = NULL;
    cal->p_free_vals = NULL;
    cal->n_elements = (size_t)0;

    /*  And no error has occurred, yet.                                       */
    cal->error_message = NULL;
    cal->error_occurred = tmpl_False;
}
/*  End of rssringoccs_CalCSV_Init.                                           */
