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
 *      Returns the y component of a three dimensional vector/spacial point.  *
 *      That is, given (x, y, z), return y.                                   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and three-vector typedef found here.                   */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for returning the y component of a three dimensional vector.     */
double rssringoccs_ThreeVector_Y(rssringoccs_ThreeVector P)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double y;

    /*  rssringoccs_ThreeVector is a struct consisting of a single double     *
     *  array dat with three entries. The y component corresponds to the      *
     *  first entry. Retrieve this and return.                                */
    y = P.dat[1];
    return y;
}
/*  End of rssringoccs_ThreeVector_Y.                                         */
