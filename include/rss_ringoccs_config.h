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

#define __RSS_RINGOCCS_USE_LOG_ALGORITHM__ 0

/*  If long double precision (80, 96, or 128 bit format, 24 decimal accuracy) *
 *  is needed for calculations but your platform lacks sinl or is strictly    *
 *  C89 compliant, rss_ringoccs provides simple yet accurate code for         *
 *  computing the various trig functions. When tested against sinl on a       *
 *  massive array of 10^8 points evenly distributed throughout -10 to 10 we   *
 *  get the following times (in seconds):                                     *
 *      C99:          4.586321                                                *
 *      rss_ringoccs: 5.659821                                                *
 *  Both outputs had long double precision and the maximum absolute           *
 *  difference was ~ 10^-16. By default we assume you want the standard       *
 *  library functions. If you do not have sinl, rss_ringoccs creates a sinl   *
 *  function via:                                                             *
 *      sinl(x) = (long double)sin((double)x)                                 *
 *  That is, we simply cast the input and use the standard sin function. This *
 *  will lose precision, but is faster. If you would like to implement our    *
 *  algorithms, set __RSS_RINGOCCS_USE_TRIG_ALGORITHMS__ to 1. The algorithms *
 *  use lookup tables in combination with Taylor series and a few trig        *
 *  identities. sin(x) is computed as follows:                                *
 *                                                                            *
 *      sin(x) = sin(y + 2pi k)                                               *
 *             = sin(y)                                                       *
 *             = sin(y' + dy)                                                 *
 *             = sin(y')cos(dy) + cos(y')sin(dy)                              *
 *                                                                            *
 *  So we get x into the range [-pi, pi] and write y = y'+dy where            *
 *  y' = 0.01*floor(100*y) and dy = y-y'. So the y' take on a small amount of *
 *  possible values and cos(y') and sin(y') are computed via lookup tables.   *
 *  dy is small and can be accurately computed via a Taylor series using very *
 *  few terms. Cosine is defined similarly.                                   */
#define __RSS_RINGOCCS_USE_TRIG_ALGORITHMS__ 0

#endif
/*  End of include guard.                                                     */
