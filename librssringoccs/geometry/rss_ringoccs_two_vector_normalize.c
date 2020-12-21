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
 *      Normalize a non-zero vector to have length 1.                         *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  NaN is defined here.                                                      */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Function prototype and two-vector typedef found here.                     */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for normalizing a non-zero vector to length 1.                   */
rssringoccs_TwoVector
rssringoccs_TwoVector_Normalize(rssringoccs_TwoVector P)
{
    double norm, rcpr_norm, x, y, x_hat, y_hat;
    rssringoccs_TwoVector P_normalized;

    norm = rssringoccs_TwoVector_Euclidean_Norm(P);

    if (norm == 0.0)
    {
        x_hat = rssringoccs_NaN;
        y_hat = rssringoccs_NaN;
    }
    else
    {
        x = rssringoccs_TwoVector_X(P);
        y = rssringoccs_TwoVector_Y(P);
        rcpr_norm = 1.0/norm;
        x_hat = x*rcpr_norm;
        y_hat = y*rcpr_norm;
    }
    /*  End of if (norm == 0.0).                                              */

    P_normalized = rssringoccs_TwoVector_Rect(x_hat, y_hat);
    return P_normalized;
}
/*  End of rssringoccs_TwoVector_Normalize.                                   */
