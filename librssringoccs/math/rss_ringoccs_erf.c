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
 *                            rss_ringoccs_erf                                *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the error function Erf.                  *
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
 *  Date:       December 10, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/10 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file where the prototypes for these functions are defined.         */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Erf is not required in C89, so we provide the algorithm for               *
 *  double, float, and long double inputs.                                    */

/*  Single precision error function (erff equivalent).                        */
float rssringoccs_Float_Erf(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float erf, erfc;

    /*  The error function can be computed by the complimentary error         *
     *  Erfc via Erf(x) = 1 - Erfc(x). We'll use this.                        */
    erfc = rssringoccs_Float_Erfc(x);
    erf = 1.0F - erfc;
    return erf;
}
/*  End of rssringoccs_Float_Erf.                                             */

/*  Double precision error function (erf equivalent).                         */
double rssringoccs_Double_Erf(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double erf, erfc;

    /*  The error function can be computed by the complimentary error         *
     *  Erfc via Erf(x) = 1 - Erfc(x). We'll use this.                        */
    erfc = rssringoccs_Double_Erfc(x);
    erf = 1.0 - erfc;
    return erf;
}
/*  End of rssringoccs_Double_Erf.                                            */


/*  Long double precision error function (erfl equivalent).                   */
long double rssringoccs_LDouble_Erf(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double erf, erfc;

    /*  The error function can be computed by the complimentary error         *
     *  Erfc via Erf(x) = 1 - Erfc(x). We'll use this.                        */
    erfc = rssringoccs_LDouble_Erfc(x);
    erf = 1.0L - erfc;
    return erf;
}
/*  End of rssringoccs_Double_Erf.                                            */
