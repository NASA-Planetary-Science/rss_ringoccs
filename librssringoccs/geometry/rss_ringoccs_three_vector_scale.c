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
 *      Function for multiplying a three vector by a real number.             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and three-vector typedef found here.                   */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Multiply a three vector by a real number.                                 */
rssringoccs_ThreeVector
rssringoccs_ThreeVector_Scale(double a, rssringoccs_ThreeVector P)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double Px, Py, Pz, x, y, z;
    rssringoccs_ThreeVector scale;

    /*  Extract the x, y, and z components from P.                            */
    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    /*  Scalar multiplication is done component-wise, so compute this.        */
    x = a*Px;
    y = a*Py;
    z = a*Pz;

    /*  Use rssringoccs_ThreeVector_Rect to create the output and return.     */
    scale = rssringoccs_ThreeVector_Rect(x, y, z);
    return scale;
}
/*  End of rssringoccs_ThreeVector_Scale.                                     */
