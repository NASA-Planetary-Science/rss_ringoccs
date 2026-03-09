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
 *      Provides tools for creating history objects.                          *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       June 14, 2021                                                 *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef RSS_RINGOCCS_HISTORY_H
#define RSS_RINGOCCS_HISTORY_H

/*  rssringoccs_History struct typedef provided here.                         */
#include <rss_ringoccs/include/types/rss_ringoccs_history.h>

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Date_To_Rev                                               *
 *  Purpose:                                                                  *
 *      Converts year + day number to Cassini Rev number.                     *
 *  Arguments:                                                                *
 *      year (const unsigned int):                                            *
 *          The year for the rev.                                             *
 *      doy (const unsigned int):                                             *
 *          The day of year (or day number).                                  *
 *  Outputs:                                                                  *
 *      rev (const char *):                                                   *
 *          The rev number as a string, like "007RI".                         *
 ******************************************************************************/
extern const char *
rssringoccs_Date_To_Rev(const unsigned int year, const unsigned int doy);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_History_Print                                             *
 *  Purpose:                                                                  *
 *      Prints a given history to stdout.                                     *
 *  Arguments:                                                                *
 *      history (const rssringoccs_History * const):                          *
 *          The history that is to be printed.                                *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void
rssringoccs_History_Print(const rssringoccs_History * const history);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_History_Init                                              *
 *  Purpose:                                                                  *
 *      Initialize a history object. All members except the input_vars and    *
 *      input_kwds variables are set to their correct values for the history. *
 *  Arguments:                                                                *
 *      history (rssringoccs_History * const):                                *
 *          The history that is to be initialized.                            *
 *  Outputs:                                                                  *
 *      None (void).                                                          *
 ******************************************************************************/
extern void rssringoccs_History_Init(rssringoccs_History * const history);

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Version                                                   *
 *  Purpose:                                                                  *
 *      Returns the current version of rss_ringoccs as a string.              *
 *  Arguments:                                                                *
 *      None (void).                                                          *
 *  Outputs:                                                                  *
 *      version (const char *):                                               *
 *          The version number as a string (i.e. "1.3").                      *
 ******************************************************************************/
extern const char *rssringoccs_Version(void);

#endif
/*  End of include guard.                                                     */
