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
 *      Given two real numbers r and theta, return the vector (x, y) where    *
 *      x = r cos(theta) and y = r sin(theta). That is, return the polar      *
 *      representation of (r, theta) in the plane.                            *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Trig functions found here.                                                */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Function prototype and two-vector typedef found here.                     */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for returning the the polar representation of (r, theta).        */
rssringoccs_TwoVector rssringoccs_TwoVector_Polar(double r, double theta)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_TwoVector P;
    double x, y;

    /*  Compute the Cartesian coordinates of (r, theta) using cosine and sine.*/
    x = r * rssringoccs_Double_Cos(theta);
    y = r * rssringoccs_Double_Sin(theta);

    /*  Set the zeroth entry of P.dat to x and the first entry to y.          */
    P.dat[0] = x;
    P.dat[1] = y;
    return P;
}
/*  End of rssringoccs_TwoVector_Polar.                                       */
