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
 *      Computes the cross product of two vectors.                            *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and three-vector typedef found here.                   */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for computing the cross product of two vectors.                  */
rssringoccs_ThreeVector
rssringoccs_Cross_Product(rssringoccs_ThreeVector P, rssringoccs_ThreeVector Q)
{
    /*  Declare necessary variables.                                          */
    double Px, Py, Pz, Qx, Qy, Qz, x, y, z;
    rssringoccs_ThreeVector cross;

    /*  Extract the x, y, and z components from P.                            */
    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    /*  Extract the x, y, and z components from Q.                            */
    Qx = rssringoccs_ThreeVector_X(Q);
    Qy = rssringoccs_ThreeVector_Y(Q);
    Qz = rssringoccs_ThreeVector_Z(Q);

    /*  Compute the components of the cross product PxQ.                      */
    x = Py*Qz - Pz*Qy;
    y = Pz*Qx - Px*Qz;
    z = Px*Qy - Py*Qx;

    /*  Use rssringoccs_ThreeVector_Rect to create the output and return.     */
    cross = rssringoccs_ThreeVector_Rect(x, y, z);
    return cross;
}
/*  End of rssringoccs_Cross_Product.                                         */
