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
 *                        rss_ringoccs_complex_dist                           *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the function f(z,w) = |z-w|.             *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CFloat_Dist:                                              *
 *      rssringoccs_CDouble_Dist:                                             *
 *      rssringoccs_CLDouble_Dist:                                            *
 *  Purpose:                                                                  *
 *      Computes the distance between two complex numbers:                    *
 *                                                                            *
 *          dist(z, w) = dist(a + ib, c + id)                                 *
 *                     = sqrt((c-a)^2 + (d-b)^2)                              *
 *  Arguments:                                                                *
 *      z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          A complex number.                                                 *
 *      w (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          Another complex number.                                           *
 *  Output:                                                                   *
 *      dist (float/double/long double):                                      *
 *          The distance between z and w.                                     *
 *  Method:                                                                   *
 *      Treat the points as elements of the Euclidean plane and use           *
 *      the Pythagorean formula.                                              *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 *  2.) rss_ringoccs_math.h:                                                  *
 *          Header containing various math functions like the square root.    *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code is compatible with the C89/C90 standard. The setup script that  *
 *  is used to compile this in config_librssringoccs.sh uses gcc and has the  *
 *  -pedantic and -std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 23, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/23 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file containing the square root functions.                         */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Complex routines and data types defined here.                             */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  Single precision distance function for complex variables.                 */
float rssringoccs_CFloat_Dist(rssringoccs_ComplexFloat z0,
                              rssringoccs_ComplexFloat z1)
{
    /*  Declare necessary variables C89 requires this at the top.             */
    float x0, y0, x1, y1, dx, dy, dist;

    /*  Extract the real and imaginary parts from z0.                         */
    x0 = rssringoccs_CFloat_Real_Part(z0);
    y0 = rssringoccs_CFloat_Imag_Part(z0);

    /*  Extract the real and imaginary parts from z1.                         */
    x1 = rssringoccs_CFloat_Real_Part(z1);
    y1 = rssringoccs_CFloat_Imag_Part(z1);

    /*  Compute the x and y differences from z0 and z1.                       */
    dx = x1-x0;
    dy = y1-y0;

    /*  Use the Pythagorean formula to compute the distance and return.       */
    dist = rssringoccs_Float_Sqrt(dx*dx + dy*dy);
    return dist;
}
/*  End of rssringoccs_CFloat_Dist.                                           */

/*  Double precision distance function for complex variables.                 */
double rssringoccs_CDouble_Dist(rssringoccs_ComplexDouble z0,
                                rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables C89 requires this at the top.             */
    double x0, y0, x1, y1, dx, dy, dist;

    /*  Extract the real and imaginary parts from z0.                         */
    x0 = rssringoccs_CDouble_Real_Part(z0);
    y0 = rssringoccs_CDouble_Imag_Part(z0);

    /*  Extract the real and imaginary parts from z1.                         */
    x1 = rssringoccs_CDouble_Real_Part(z1);
    y1 = rssringoccs_CDouble_Imag_Part(z1);

    /*  Compute the x and y differences from z0 and z1.                       */
    dx = x1-x0;
    dy = y1-y0;

    /*  Use the Pythagorean formula to compute the distance and return.       */
    dist = rssringoccs_Double_Sqrt(dx*dx + dy*dy);
    return dist;
}
/*  End of rssringoccs_CDouble_Dist.                                          */

/*  Long double precision distance function for complex variables.            */
long double rssringoccs_CLDouble_Dist(rssringoccs_ComplexLongDouble z0,
                                      rssringoccs_ComplexLongDouble z1)
{
    /*  Declare necessary variables C89 requires this at the top.             */
    long double x0, y0, x1, y1, dx, dy, dist;

    /*  Extract the real and imaginary parts from z0.                         */
    x0 = rssringoccs_CLDouble_Real_Part(z0);
    y0 = rssringoccs_CLDouble_Imag_Part(z0);

    /*  Extract the real and imaginary parts from z1.                         */
    x1 = rssringoccs_CLDouble_Real_Part(z1);
    y1 = rssringoccs_CLDouble_Imag_Part(z1);

    /*  Compute the x and y differences from z0 and z1.                       */
    dx = x1-x0;
    dy = y1-y0;

    /*  Use the Pythagorean formula to compute the distance and return.       */
    dist = rssringoccs_LDouble_Sqrt(dx*dx + dy*dy);
    return dist;
}
/*  End of rssringoccs_CLDouble_Dist.                                         */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here then we have C99 complex.h support and can use the         *
 *  arithmetic combined with the cabs function.                               */

/*  Single precision distance function for complex variables.                 */
float rssringoccs_CFloat_Dist(rssringoccs_ComplexFloat z0,
                               rssringoccs_ComplexFloat z1)
{
    return cabsf(z0 - z1);
}
/*  End of rssringoccs_CFloat_Dist.                                           */

/*  Double precision distance function for complex variables.                 */
double rssringoccs_CDouble_Dist(rssringoccs_ComplexDouble z0,
                                rssringoccs_ComplexDouble z1)
{
    return cabs(z0 - z1);
}
/*  End of rssringoccs_CDouble_Dist.                                          */

/*  Long double precision distance function for complex variables.            */
long double rssringoccs_CLDouble_Dist(rssringoccs_ComplexLongDouble z0,
                                      rssringoccs_ComplexLongDouble z1)
{
    return cabsl(z0 - z1);
}
/*  End of rssringoccs_CLDouble_Dist.                                         */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */
