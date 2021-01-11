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
 *      Computes the Euclidean norm of a three dimensional vector using the   *
 *      Pythagorean formula:                                                  *
 *          ||(x, y, z)|| = sqrt(x^2 + y^2 + z^2)                             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 21, 2020                                             *
 ******************************************************************************/

/*  Square root function found here.                                          */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Function prototype and three-vector typedef found here.                   */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Function for computing the length of three dimensional vectors.           */
double rssringoccs_ThreeVector_Euclidean_Norm(rssringoccs_ThreeVector P)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double x, y, z, norm;

    /*  Extract the x, y, and z components from P.                            */
    x = rssringoccs_ThreeVector_X(P);
    y = rssringoccs_ThreeVector_Y(P);
    z = rssringoccs_ThreeVector_Z(P);

    /*  Use the Pythagorean formula to compute the norm and return.           */
    norm = rssringoccs_Double_Sqrt(x*x + y*y + z*z);
    return norm;
}
/*  End of rssringoccs_Euclidean_Norm_3D.                                     */
