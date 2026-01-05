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
 *      Provides a struct for the data in a CAL.TAB file.                     *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 5, 2026                                               *
 ******************************************************************************/

/*  Include guard to prevent including this file.                             */
#ifndef RSS_RINGOCCS_TYPES_CALCSV_H
#define RSS_RINGOCCS_TYPES_CALCSV_H

/*  Booleans found here.                                                      */
#include <libtmpl/include/tmpl_bool.h>

/*  History object typedef is here. Each CSV object gets its own history.     */
#include <rss_ringoccs/include/types/rss_ringoccs_history.h>

/*  size_t typedef provided here.                                             */
#include <stddef.h>

/*  Data structure for the CAL.TAB files on the PDS.                          */
typedef struct rssringoccs_CalCSV_Def {
    double *t_oet_spm_vals;
    double *f_sky_pred_vals;
    double *f_sky_resid_fit_vals;
    double *p_free_vals;
    size_t n_elements;
    rssringoccs_History *history;
    tmpl_Bool error_occurred;
    const char *error_message;
} rssringoccs_CalCSV;

#endif
/*  End of include guard.                                                     */
