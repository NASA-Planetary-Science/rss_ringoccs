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
 *                             rss_ringoccs_abs                               *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the absolute value.                      *
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

/*  Since the absolute value function is very simple, we simply provide the   *
 *  algorithms here rather than pass the arguments to fabs, fabsf, or fabsfl. *
 *  There is essentially no time difference. Using gcc with -O2 optimization  *
 *  on an array of 1 million random elements in [-1, 1] gave the following    *
 *  times (in seconds):                                                       *
 *      fabs:         0.003328                                                *
 *      rss_ringoccs: 0.003743                                                *
 *  -O3 optimization gives:                                                   *
 *      fabs:         0.003409                                                *
 *      rss_ringoccs: 0.003493                                                *
 *  So, no real difference. These times were computed on an iMac 2017 3.4GHz  *
 *  quad-core running MacOS Catalina 10.15.7. Converting a long double to a   *
 *  double will lose precision, hence the reason we provide this simple code. */

/*  Single precision absolute value function (fabsf equivalent).              */
float rssringoccs_Float_Abs(float x)
{
    /*  If x is positive return it, otherwise return its negative.            */
    if (x >= 0.0F)
        return x;
    else
        return -x;
}
/*  End of rssringoccs_Float_Abs.                                             */

/*  Double precision absolute value function (fabs equivalent).               */
double rssringoccs_Double_Abs(double x)
{
    /*  If x is positive return it, otherwise return its negative.            */
    if (x >= 0.0)
        return x;
    else
        return -x;
}
/*  End of rssringoccs_Double_Abs.                                            */

/*  Long double precision absolute value function (fabsl equivalent).         */
long double rssringoccs_LDouble_Abs(long double x)
{
    /*  If x is positive return it, otherwise return its negative.            */
    if (x >= 0.0L)
        return x;
    else
        return -x;
}
/*  End of rssringoccs_LDouble_Abs.                                           */
