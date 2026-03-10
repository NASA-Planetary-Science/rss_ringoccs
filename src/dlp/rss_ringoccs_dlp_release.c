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

/*  Header file with the DLPObj typedef and helper routines.                  */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  Function for safely decrementing the reference counter for a DLP object.  */
void rssringoccs_DLP_Release(rssringoccs_DLPObj *dlp)
{
    /*  If the input is NULL, there is no counter to decrement. Return.       */
    if (!dlp)
        return;

    /*  If the reference counter is already zero, then this DLP holds no data.*
     *  We do not need to modify the reference count or free anything. Return.*/
    if (dlp->reference_count == 0)
        return;

    /*  Otherwise decrement the count.                                        */
    --dlp->reference_count;

    /*  If the count is now zero, then nothing is referencing this data       *
     *  anymore and it needs to be free'd. Check for this.                    */
    if (dlp->reference_count == 0)
        rssringoccs_DLP_Destroy(&dlp);
}
/*  End of rssringoccs_DLP_Release.                                           */
