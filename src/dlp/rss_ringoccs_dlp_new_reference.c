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
 *  Date:       March 3, 2026                                                 *
 ******************************************************************************/

/*  Header file with the DLPObj typedef.                                      */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/*  NULL provided here.                                                       */
#include <stddef.h>

/*  Function prototype / forward declaration.                                 */
extern rssringoccs_DLPObj *
rssringoccs_DLP_New_Reference(rssringoccs_DLPObj * const dlp);

/*  Function for safely incrementing the reference counter for a DLP object.  */
rssringoccs_DLPObj *
rssringoccs_DLP_New_Reference(rssringoccs_DLPObj * const dlp)
{
    /*  If the input is NULL, there is no DLP to reference. Return NULL.      */
    if (!dlp)
        return NULL;

    /*  A DLP object that encountered an error should not be referenced. Do   *
     *  not increment the counter, instead return NULL.                       */
    if (dlp->error_occurred)
        return NULL;

    /*  Otherwise we have a valid DLP object to work with. To avoid free'ing  *
     *  the data while some instance of it is needed, increment the reference *
     *  counter.                                                              */
    ++dlp->reference_count;

    /*  The new reference is simply the original pointer with the reference   *
     *  count increased. It is the callers responsibility to "release" this   *
     *  reference using rssringoccs_DLP_Release.                              */
    return dlp;
}
/*  End of rssringoccs_DLP_New_Reference.                                     */
