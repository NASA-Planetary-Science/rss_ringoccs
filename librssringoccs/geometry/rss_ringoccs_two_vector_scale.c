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
 *      Source code for multiplying a two dimensional vector by a real number.*
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Function prototype and two-vector typedef found here.                     */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for multipyling a two-vector by a real number.                   */
rssringoccs_TwoVector
rssringoccs_TwoVector_Scale(double a, rssringoccs_TwoVector P)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double Px, Py, x, y;
    rssringoccs_TwoVector scale;

    /*  Extract the x and y components from P.                                */
    Px = rssringoccs_TwoVector_X(P);
    Py = rssringoccs_TwoVector_Y(P);

    /*  Multiplying a vector by a scalar simply multiplies the entries        *
     *  component-wise. Compute this.                                         */
    x = a*Px;
    y = a*Py;

    /*  Use rssringoccs_TwoVector_Rect to create the output and return.       */
    scale = rssringoccs_TwoVector_Rect(x, y);
    return scale;
}
/*  End of rssringoccs_TwoVector_Scale.                                       */
