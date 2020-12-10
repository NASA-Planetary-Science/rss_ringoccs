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
 *                          rss_ringoccs_copysign                             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the copysign function defined in C99.    *
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
 *  Date:       December 8, 2020                                              *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/10 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file where the prototypes for these functions are defined.         */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Copysign is not required in C89, so we provide the algorithm for          *
 *  double, float, and long double inputs.                                    */

/*  Float precision copysign function.                                        */
float rssringoccs_Float_Copysign(float x, float y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float out;

    /*  If y is negative, compute -|x|.                                       */
    if (y < 0.0F)
        out = -rssringoccs_Float_Abs(x);

    /*  If y is positive, compute |x|.                                        */
    else if (0.0F < y)
        out = rssringoccs_Float_Abs(x);

    /*  And lastly, if y is zero, return zero.                                */
    else
        out = 0.0F;

    return out;
}
/*  End of rssringoccs_Float_Copysign.                                        */

/*  Double precision copysign function.                                       */
double rssringoccs_Double_Copysign(double x, double y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double out;

    /*  If y is negative, compute -|x|.                                       */
    if (y < 0.0)
        out = -rssringoccs_Double_Abs(x);

    /*  If y is positive, compute |x|.                                        */
    else if (0.0 < y)
        out = rssringoccs_Double_Abs(x);

    /*  And lastly, if y is zero, return zero.                                */
    else
        out = 0.0;

    return out;
}
/*  End of rssringoccs_Double_Copysign.                                       */

/*  Long double precision copysign function.                                  */
long double rssringoccs_LDouble_Copysign(long double x, long double y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double out;

    /*  If y is negative, compute -|x|.                                       */
    if (y < 0.0L)
        out = -rssringoccs_LDouble_Abs(x);

    /*  If y is positive, compute |x|.                                        */
    else if (0.0L < y)
        out = rssringoccs_LDouble_Abs(x);

    /*  And lastly, if y is zero, return zero.                                */
    else
        out = 0.0L;

    return out;
}
/*  End of rssringoccs_LDouble_Copysign.                                      */
