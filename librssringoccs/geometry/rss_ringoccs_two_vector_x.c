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
 *  Purpose:                                                                  *
 *      Returns the x component of a two dimensional vector/planar point.     *
 *      That is, given (x, y), return x.                                      *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and two-vector typedef found here.                     */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for returning the x component of a two dimensional vector.       */
double rssringoccs_TwoVector_X(rssringoccs_TwoVector P)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double x;

    /*  rssringoccs_TwoVector is defined the same way                         *
     *  rssringoccs_ComplexDouble is. It is a struct containing a single      *
     *  double array dat consisting of two entries. The zeroth entry          *
     *  corresponds to the x component. Retrieve this and return.             */
    x = P.dat[0];
    return x;
}
/*  End of rssringoccs_TwoVector_X.                                           */
