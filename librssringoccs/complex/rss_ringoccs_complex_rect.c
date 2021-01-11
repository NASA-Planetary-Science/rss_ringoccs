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
 *                        rss_ringoccs_complex_rect                           *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for creating complex numbers.                *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 1, 2020                                              *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  Create single precision complex numbers in Cartesian coordinates.         */
rssringoccs_ComplexFloat rssringoccs_CFloat_Rect(float x, float y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat z;

    /*  Simply set the array dat inside z to {x, y} and return.               */
    z.dat[0] = x;
    z.dat[1] = y;
    return z;
}
/*  End of rssringoccs_CFloat_Rect.                                     */

/*  Create double precision complex numbers in Cartesian coordinates.         */
rssringoccs_ComplexDouble rssringoccs_CDouble_Rect(double x, double y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z;

    /*  Simply set the array dat inside z to {x, y} and return.               */
    z.dat[0] = x;
    z.dat[1] = y;
    return z;
}
/*  End of rssringoccs_CDouble_Rect.                                    */

/*  Create long double precision complex numbers in Cartesian coordinates.    */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Rect(long double x, long double y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble z;

    /*  Simply set the array dat inside z to {x, y} and return.               */
    z.dat[0] = x;
    z.dat[1] = y;
    return z;
}
/*  End of rssringoccs_CLDouble_Rect.                                */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  Create single precision complex numbers in Cartesian coordinates.         */
rssringoccs_ComplexFloat rssringoccs_CFloat_Rect(float x, float y)
{
    return x + _Complex_I*y;
}
/*  End of rssringoccs_CFloat_Rect.                                     */

/*  Create double precision complex numbers in Cartesian coordinates.         */
rssringoccs_ComplexDouble rssringoccs_CDouble_Rect(double x, double y)
{
    return x + _Complex_I*y;
}
/*  End of rssringoccs_CDouble_Rect.                                    */

/*  Create long double precision complex numbers in Cartesian coordinates.    */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Rect(long double x, long double y)
{
    return x + _Complex_I*y;
}
/*  End of rssringoccs_CLDouble_Rect.                                */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */
