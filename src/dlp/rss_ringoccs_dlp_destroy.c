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

void rssringoccs_DLP_Destroy(rssringoccs_DLPObj ** const dlp)
{
    if (!dlp)
        return;

    if (!(*dlp))
        return;

    rssringoccs_DLP_Destroy_Members(*dlp);
    TMPL_FREE(*dlp);
}
/*  End of rssringoccs_DLP_Destroy.                                           */
