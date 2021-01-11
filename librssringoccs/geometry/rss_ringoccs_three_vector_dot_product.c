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
 *      Compute the Euclidean dot product (a,b,c) . (x,y,z) = ax + by + cz.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and thre-vector typedef found here.                    */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for computing the dot product of 2 three-vectors.                */
double rssringoccs_ThreeVector_Dot_Product(rssringoccs_ThreeVector P,
                                           rssringoccs_ThreeVector Q)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double Px, Py, Pz, Qx, Qy, Qz, dot_prod;

    /*  Extract the x, y, and z components from P.                            */
    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    /*  Extract the x, y, and z components from Q.                            */
    Qx = rssringoccs_ThreeVector_X(Q);
    Qy = rssringoccs_ThreeVector_Y(Q);
    Qz = rssringoccs_ThreeVector_Z(Q);

    /*  Use the Euclidean dot product formula and return.                     */
    dot_prod = Px*Qx + Py*Qy + Pz*Qz;
    return dot_prod;
}
/*  End of rssringoccs_ThreeVector_Dot_Product.                               */
