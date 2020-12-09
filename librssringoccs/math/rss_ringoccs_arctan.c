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
 *                           rss_ringoccs_arctan                              *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the inverse tangent function.            *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 1, 2020                                              *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/08 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file where the prototypes for these functions are defined.         */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  The "double" version of atan is defined in both C89 and C99 math.h so we  *
 *  only need to alias this function.                                         */

/*  Double precision arctan function (atan equivalent).                       */
double rssringoccs_Double_Arctan(double x)
{
    return atan(x);
}
/*  End of rssringoccs_Double_Arctan.                                         */

/*  Check if we have C99 math.h or not.                                       */
#if __HAS_C99_MATH_H__ == 0

/*  C89 math.h does not have atanf or atanl, so we'll need to provide these   *
 *  to make the code forward compatible.                                      */

/*  Single precision arctan function (atanf equivalent).                      */
float rssringoccs_Float_Arctan(float x)
{
    /*  Cast "x" as a double to atan, and cast the output as a float.         */
    return (float)atan((double)x);
}
/*  End of rssringoccs_Float_Arctan.                                          */

/*  Long double precision arctan function (atanl equivalent).                 */
long double rssringoccs_LDouble_Arctan(long double x)
{
    /*  Cast "x" as a double to atan, and cast the output as a long double.   */
    return (long double)atan((double)x);
}
/*  End of rssringoccs_LDouble_Arctan.                                        */

#else
/*  Else statement for #if __HAS_C99_MATH_H__ == 0.                           */

/*  C99 provides float and long double support for their math functions, so   *
 *  simply use to these.                                                      */

/*  Single precision arctan function (atanf equivalent).                      */
float rssringoccs_Float_Arctan(float x)
{
    return atanf(x);
}
/*  End of rssringoccs_Float_Arctan.                                          */

/*  Long double precision arctan function (atanl equivalent).                 */
long double rssringoccs_LDouble_Arctan(long double x)
{
    return atanl(x);
}
/*  End of rssringoccs_LDouble_Arctan.                                        */

#endif
/*  End of #if __HAS_C99_MATH_H__ == 0                                        */
