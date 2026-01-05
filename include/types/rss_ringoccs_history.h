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
 *      Provides a struct for tracking the history of a processing pipeline.  *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       January 5, 2026                                               *
 ******************************************************************************/

/*  Include guard to prevent including this file.                             */
#ifndef RSS_RINGOCCS_TYPES_HISTORY_H
#define RSS_RINGOCCS_TYPES_HISTORY_H

/*  Struct with information about how a data set was processed.               */
typedef struct rssringoccs_History_Def {
    const char *rss_ringoccs_version;
    const char *libtmpl_version;
    const char *c_version;
    const char *user_name;
    const char *host_name;
    const char *run_date;
    const char *operating_system;

    /*  To avoid having to malloc this array, it is statically sized. The     *
     *  strings come in pairs "parameter_name: the_parameter", and the array  *
     *  ends in a NULL terminator. So you can have (17 - 1) / 2 = 8 different *
     *  arguments, and similarly 8 different keywords. For the rss_ringoccs   *
     *  structs that have history data, this is always sufficient.            */
    const char *input_vars[17];
    const char *input_kwds[17];
} rssringoccs_History;

#endif
/*  End of include guard.                                                     */
