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
 *                   rss_ringoccs_stereographic_projection                    *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains the source code for the stereographic projection   *
 *      of a point on a sphere to the plane.                                  *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_geometry.h:                                              *
 *          rssringoccs_TwoVector and rssringoccs_ThreeVector are defined     *
 *          here, as well as the prototypes for the functions.                *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 23, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/11/24 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  rssringoccs_NaN is defined here.                                          */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/* rssringoccs_TwoVector and rssringoccs_ThreeVector are defined here.        */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Compute the stereographic projection of a point on a sphere. The radius   *
 *  of the sphere is given by the length of the input vector P.               */
rssringoccs_TwoVector
rssringoccs_Stereographic_Projection(rssringoccs_ThreeVector P)
{
    /*  Declare all necessary variables. C89 requires this at the top.        */
    double Px, Py, Pz, x, y, factor, norm;
    rssringoccs_TwoVector out;

    /*  Extract the X, Y, and Z components of the input vector.               */
    Px = rssringoccs_ThreeVector_X(P);
    Py = rssringoccs_ThreeVector_Y(P);
    Pz = rssringoccs_ThreeVector_Z(P);

    /*  We will perform stereographic projection on the sphere centered at    *
     *  the origin which contains the given point. To do this we'll first     *
     *  need to compute the Euclidean norm of this point.                     */
    norm = rssringoccs_ThreeVector_Euclidean_Norm(P);

    /*  In the special case of P = 0, stereographic projection is not defined *
     *  and we can't say it's the "point at infinity." We'll just return      *
     *  (0, 0) in this instance.                                              */
    if (norm == 0.0)
    {
        x = 0.0;
        y = 0.0;
    }

    /*  The next extreme case is when Pz = norm. That is, our point is the    *
     *  north pole. In this instance stereographic projection takes us to     *
     *  "the point at infinity" so we'll return (infinity, infinity).         */
    else if (Pz == norm)
    {
        x = rssringoccs_Infinity;
        y = rssringoccs_Infinity;
    }

    /*  Otherwise we have a "normal" case and can use stereographic           *
     *  projection. This is given by the following equations.                 */
    else
    {
        factor = norm / (norm - Pz);
        x = factor * Px;
        y = factor * Py;
    }

    /*  Use rssringoccs_TwoVector_Rect to create a TwoVector from x and y     *
     *  and return.                                                           */
    out = rssringoccs_TwoVector_Rect(x, y);
    return out;
}
