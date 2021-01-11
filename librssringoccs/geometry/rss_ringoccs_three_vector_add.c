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
 *      Source code for adding three dimensional vectors.                     *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and three-vector typedef found here.                   */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for adding 2 three-dimensional vectors.                          */
rssringoccs_ThreeVector
rssringoccs_ThreeVector_Add(rssringoccs_ThreeVector P,
                            rssringoccs_ThreeVector Q)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double Px, Py, Pz, Qx, Qy, Qz, x, y, z;
    rssringoccs_ThreeVector sum;

    /*  Extract the x, y, and z components from P.                            */
    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    /*  Extract the x, y, and z components from Q.                            */
    Qx = rssringoccs_ThreeVector_X(Q);
    Qy = rssringoccs_ThreeVector_Y(Q);
    Qz = rssringoccs_ThreeVector_Z(Q);

    /*  The sum of two vectors simply adds their components together.         */
    x = Px + Qx;
    y = Py + Qy;
    z = Pz + Qz;

    /*  Use rssringoccs_ThreeVector_Rect to create the output and return.     */
    sum = rssringoccs_ThreeVector_Rect(x, y, z);
    return sum;
}
/*  End of rssringoccs_ThreeVector_Add.                                       */
