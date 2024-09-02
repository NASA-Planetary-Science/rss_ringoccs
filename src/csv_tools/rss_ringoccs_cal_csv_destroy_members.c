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

/*  free is found here, as is NULL.                                           */
#include <stdlib.h>

/*  rssringoccs_CalCSV typedef here, and function prototype given.            */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

/*  Check if this macro name has been defined. Clear it if so, as we'll be    *
 *  using this macro to free the variables in the CAL struct.                 */
#ifdef DESTROY_CAL_VAR
#undef DESTROY_CAL_VAR
#endif

/*  Macro for freeing and nullifying the members of the cal CSV structs.      */
#define DESTROY_CAL_VAR(var) if (var != NULL){free(var); var = NULL;}

/*  Free's all members of an rssringoccs_CalCSV pointer except the            *
 *  error_message. Members are set to NULL after freeing.                     */
void rssringoccs_CalCSV_Destroy_Members(rssringoccs_CalCSV *cal)
{
    /*  If the pointer is NULL, there's nothing to do. Simply return.         */
    if (!cal)
        return;

    /*  Destroy every variable except the error_message.                      */
    DESTROY_CAL_VAR(cal->t_oet_spm_vals)
    DESTROY_CAL_VAR(cal->f_sky_pred_vals)
    DESTROY_CAL_VAR(cal->f_sky_resid_fit_vals)
    DESTROY_CAL_VAR(cal->p_free_vals)
}
/*  End of rssringoccs_CalCSV_Destroy_Members.                                */

/*  Undefine everything in case someone wants to #include this file.          */
#undef DESTROY_CAL_VAR
