/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *                            rss_ringoccs_bool                               *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Define bool, false, and true for C89/C90 compilers. If your compiler  *
 *      supports stdbool (i.e., C99, C11, or C18), then the built-in _Bool    *
 *      is compatible with what we've defined here.                           *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 12, 2020                                            *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_BOOL_H__
#define __RSS_RINGOCCS_BOOL_H__

/*  The C99, C11, and C18 standards require support for booleans, but it is   *
 *  not available in C89/C90. We've typedef'd Booleans as rssringoccs_Bool to *
 *  avoid naming conflicts with C99 and higher Booleans. If you compiler      *
 *  supports Booleans, you should be able to use stdbools with rss_ringoccs.  */

/*  We prepend rssringoccs onto False, True, and Bool to avoid name conflicts.*/
typedef enum {rssringoccs_False, rssringoccs_True} rssringoccs_Bool;

#endif
/*  End of include guard.                                                     */
