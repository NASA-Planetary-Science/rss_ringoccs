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
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 5, 2020                                               *
 ******************************************************************************/
#include <libtmpl/include/compat/tmpl_free.h>
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for freeing all member of a tau object except the error message. */
void rssringoccs_Tau_Destroy_Members(rssringoccs_TAUObj * const tau)
{
    /*  If the input pointer is NULL, do not try to access it. Just return.   */
    if (!tau)
        return;

    /*  Use the TMPL_FREE macro to free and Nullify all pointers.             */
    TMPL_FREE(tau->T_in);
    TMPL_FREE(tau->T_out);
    TMPL_FREE(tau->T_fwd);
    TMPL_FREE(tau->F_km_vals);
    TMPL_FREE(tau->k_vals);
    TMPL_FREE(tau->w_km_vals);
    TMPL_FREE(tau->tau_threshold_vals);

    /*  Lastly, release the DLP object and set the pointer to NULL.           */
    rssringoccs_DLP_Release(tau->dlp);
    tau->dlp = NULL;
}
/*  End of rssringoccs_Tau_Destroy_Members.                                   */
