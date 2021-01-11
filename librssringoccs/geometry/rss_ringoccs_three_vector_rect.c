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
 *      Given three real numbers x, y, and z, return the vector (x, y, z).    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and three-vector typedef found here.                   */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for returning the point (x, y, z) given three doubles x, y, z.   */
rssringoccs_ThreeVector
rssringoccs_ThreeVector_Rect(double x, double y, double z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ThreeVector P;

    /*  Set the zeroth entry of P.dat to x, the first to y, and second to z.  */
    P.dat[0] = x;
    P.dat[1] = y;
    P.dat[2] = z;
    return P;
}
/*  End of rssringoccs_ThreeVector_Rect.                                      */
