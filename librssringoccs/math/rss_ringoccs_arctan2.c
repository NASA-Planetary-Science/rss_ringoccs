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
 *                           rss_ringoccs_arctan2                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the 2-dimensional inverse tangent        *
 *      function. This computes the angle between the point (x,y) and (1,0)   *
 *      in the Cartesian plane.                                               *
 *  NOTE:                                                                     *
 *      By convention dating back to (at least) the 1970s, Arctan2 takes the  *
 *      input as (y,x), not (x,y). i.e. the first argument is the y           *
 *      component and the second argument is the x component. This is contrary*
 *      to most 2 dimensional functions that want their inputs as (x,y).      *
 *      This is probably because we are trying to compute tan^-1(y/x) but     *
 *      need to be careful about the signs of y and x, so we write            *
 *      arctan(y,x).                                                          *
 *                                                                            *
 *      This returns a number between -pi and pi, so there is a "branch cut"  *
 *      along the negative x axis. Because of this, use of this function      *
 *      in complex routines results in actual branch cuts.                    *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Float_Arctan2:                                            *
 *      rssringoccs_Double_Arctan2:                                           *
 *      rssringoccs_LDouble_Arctan2:                                          *
 *  Purpose:                                                                  *
 *      Computes arctangent of a point in the plane.                          *
 *  Arguments:                                                                *
 *      y (float/double/long double):                                         *
 *          A real number, the y component of a point in the plane.           *
 *      x (float/double/long double):                                         *
 *          A real number, the x component of a point in the plane.           *
 *  Output:                                                                   *
 *      arctan_xy (float/double/long double):                                 *
 *          The angle made between the point (x, y) and (1, 0) in radians.    *
 *  Method:                                                                   *
 *      Passes arguments to atan2f, atan2, atan2l if available, and converts  *
 *      the arguments to doubles and passes them to atan2 if not.             *
 *  Notes:                                                                    *
 *      It is not assumed you have a C99 compliant compiler. If you do and    *
 *      would like to use those library functions, set the macro              *
 *      __RSS_RINGOCCS_USING_C99_MATH_H__ to 1 in rss_ringoccs_configs.h.     *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
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
 *  -pedantic and =std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
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

/*  The "double" version of atan2 is defined in both C89 and C99 math.h so we *
 *  only need to alias this function.                                         */

/*  Double precision 2-dimensional arctan function (atan2 equivalent).        */
double rssringoccs_Double_Arctan2(double y, double x)
{
    return atan2(y, x);
}
/*  End of rssringoccs_Double_Arctan2.                                        */

/*  Check if we have C99 math.h or not.                                       */
#if __HAS_C99_MATH_H__ == 0

/*  C89 math.h does not have atan2f or atan2l, so we'll need to provide these *
 *  to make the code forward compatible.                                      */

/*  Single precision 2-dimensional arctan function (atan2f equivalent).       */
float rssringoccs_Float_Arctan2(float y, float x)
{
    /*  Cast "x" and "y" as doubles to atan2, and cast the output as a float. */
    return (float)atan2((double)y, (double)x);
}
/*  End of rssringoccs_Float_Arctan2.                                         */

/*  Long double precision 2-dimensional arctan function (atan2l equivalent).  */
long double rssringoccs_LDouble_Arctan2(long double y, long double x)
{
    /*  Cast "x" and "y" as doubles to atan2, and cast the output as a        *
     *  long double.                                                          */
    return (long double)atan2((double)y, (double)x);
}
/*  End of rssringoccs_LDouble_Arctan2.                                       */

#else
/*  Else statement for #if __HAS_C99_MATH_H__ == 0.                           */

/*  C99 provides float and long double support for their math functions, so   *
 *  simply use to these.                                                      */

/*  Single precision 2-dimensional arctan function (atan2f equivalent).       */
float rssringoccs_Float_Arctan2(float y, float x)
{
    return atan2f(y, x);
}
/*  End of rssringoccs_Float_Arctan2.                                         */

/*  Long double precision 2-dimensional arctan function (atan2l equivalent).  */
long double rssringoccs_LDouble_Arctan2(long double y, long double x)
{
    return atan2l(y, x);
}
/*  End of rssringoccs_LDouble_Arctan2.                                       */

#endif
/*  End of #if __HAS_C99_MATH_H__ == 0                                        */
