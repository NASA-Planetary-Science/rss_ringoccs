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
 *      Compute the Euclidean dot product (a,b) . (c,d) = ac + bd.            *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and two-vector typedef found here.                     */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for returning the x component of a two dimensional vector.       */
double rssringoccs_TwoVector_Dot_Product(rssringoccs_TwoVector P,
                                         rssringoccs_TwoVector Q)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double Px, Py, Qx, Qy, dot_prod;

    /*  Extract the x and y components from P.                                */
    Px = rssringoccs_TwoVector_X(P);
    Py = rssringoccs_TwoVector_Y(P);

    /*  Extract the x and y components from Q.                                */
    Qx = rssringoccs_TwoVector_X(Q);
    Qy = rssringoccs_TwoVector_Y(Q);

    /*  Use the Euclidean dot product formula and return.                     */
    dot_prod = Px*Qx + Py*Qy;
    return dot_prod;
}
/*  End of rssringoccs_TwoVector_Dot_Product.                                 */
