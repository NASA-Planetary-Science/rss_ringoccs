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
 *  Date:       February 27, 2026                                             *
 ******************************************************************************/

/*  TMPL_FREE macro provided here, free's and nullifies a pointer.            */
#include <libtmpl/include/compat/tmpl_free.h>

/*  DLP object typedef found here and helper routines provided.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  NULL macro provided here.                                                 */
#include <stddef.h>

/*  Function for destroying all of the memory in a dlp object.                */
void rssringoccs_DLP_Destroy(rssringoccs_DLPObj ** const dlp)
{
    /*  If the input pointer is NULL, do not try to access it. Just return.   */
    if (!dlp)
        return;

    /*  Similarly, if the object pointer itself is NULL, simply return.       */
    if (!(*dlp))
        return;

    /*  Free all of the memory allocated to the arrays inside the object.     */
    rssringoccs_DLP_Destroy_Members(*dlp);

    /*  Set the error members to their zero values.                           */
    (*dlp)->error_message = NULL;
    (*dlp)->error_occurred = tmpl_False;

    /*  Lastly, free and nullify the DLP object itself.                       */
    TMPL_FREE(*dlp);
}
/*  End of rssringoccs_DLP_Destroy.                                           */
