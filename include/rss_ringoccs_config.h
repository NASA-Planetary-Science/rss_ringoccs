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
 *                            rss_ringoccs_config                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Provide macros describing how to build rss_ringoccs. This includes    *
 *      the macros for using C99 complex and C99 math.h. The default build    *
 *      assumes C89 compliance and does not use any of the features from C99. *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 2, 2020                                              *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_CONFIG_H__
#define __RSS_RINGOCCS_CONFIG_H__

/*  Set this value to either 0 or 1. Setting it to zero means we'll use the   *
 *  complex data types defined in rss_ringoccs_complex.h instead of the       *
 *  ones provided by C99, setting it to 1 uses C99 complex.h. By default      *
 *  rss_ringoccs does not assume you have a C99 capable compiler and          *
 *  builds using C89/C90 (also known as ANSI C) compliant code.               */
#define __RSS_RINGOCCS_USING_COMPLEX_H__ 0

/*  There are subtle differences between the C89/C90 and C99 math.h header    *
 *  files. The C99 math.h is a superset of the C89 one. rss_ringoccs provides *
 *  the functions that are missing in rss_ringoccs_math.h. If you have a      *
 *  compiler that supports C99 math.h, then setting this macro to 1 will make *
 *  rss_ringoccs alias those functions rather than providing new algorithms.  *
 *  By default we do not assume C99 compliance.                               */
#define __RSS_RINGOCCS_USING_C99_MATH_H__ 0

#endif
/*  End of include guard.                                                     */
