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

/*  Function prototype and three-vector typedef found here.                   */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function that normalizes non-zero three dimensional vectors.              */
rssringoccs_ThreeVector
rssringoccs_ThreeVector_Normalize(rssringoccs_ThreeVector P)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double norm, rcpr_norm, x, y, z, x_hat, y_hat, z_hat;
    rssringoccs_ThreeVector P_normalized;

    /*  Get the norm of the input vector P.                                   */
    norm = rssringoccs_ThreeVector_Euclidean_Norm(P);

    /*  If the norm is zero we cannot normalize. Return NaN in this case.     */
    if (norm == 0.0)
    {
        x_hat = rssringoccs_NaN;
        y_hat = rssringoccs_NaN;
        z_hat = rssringoccs_NaN;
    }
    else
    {
        /*  Extract the x, y, and z components from P.                        */
        x = rssringoccs_ThreeVector_X(P);
        y = rssringoccs_ThreeVector_Y(P);
        z = rssringoccs_ThreeVector_Z(P);

        /*  Compute the reciprocal of the norm. Precomputing a division and   *
         *  using multiplication later is faster than repeated division.      */
        rcpr_norm = 1.0/norm;

        /*  Compute the components of the normalized vector.                  */
        x_hat = x*rcpr_norm;
        y_hat = y*rcpr_norm;
        z_hat = z*rcpr_norm;
    }
    /*  End of if (norm == 0.0).                                              */

    P_normalized = rssringoccs_ThreeVector_Rect(x_hat, y_hat, z_hat);
    return P_normalized;
}
/*  End of rssringoccs_ThreeVector_Normalize.                                 */
