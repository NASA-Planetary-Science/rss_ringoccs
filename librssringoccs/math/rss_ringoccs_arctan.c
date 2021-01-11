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
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Float_Arctan:                                             *
 *      rssringoccs_Double_Arctan:                                            *
 *      rssringoccs_LDouble_Arctan:                                           *
 *  Purpose:                                                                  *
 *      Computes arctangent of a real number:                                 *
 *          arctan(x) = tan^-1(x)                                             *
 *  Arguments:                                                                *
 *      x (float/double/long double):                                         *
 *          A real number, the argument for arctan.                           *
 *  Output:                                                                   *
 *      arctan_x (float/double/long double):                                  *
 *          The arctangent of x.                                              *
 *  Method:                                                                   *
 *      Passes arguments to atanf, atan, atanl if available, and converts the *
 *      arguments to doubles and passes them to atan if not.                  *
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
