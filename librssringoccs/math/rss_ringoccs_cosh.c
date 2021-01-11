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
 *                            rss_ringoccs_cosh                               *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the hyperbolic cosine function cosh.     *
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

/*  Cosh is not required in C89, so we provide the algorithm for              *
 *  double, float, and long double inputs.                                    */

/*  Single precision hyperbolic cosine (coshf equivalent).                    */
float rssringoccs_Float_Cosh(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float cosh_x, exp_x, exp_minus_x;

    /*  The definition of cosh(x) is [exp(x) + exp(-x)]/2, so return this. It *
     *  is computationally faster to compute exp(x) and then exp(-x) via the  *
     *  formula exp(-x) = 1/exp(x). This saves us from computing two          *
     *  exponentials at the cost of an extra division.                        */
    exp_x = rssringoccs_Float_Exp(x);
    exp_minus_x = 1.0F/exp_x;

    /*  Compute cosh from the two exponentials and return.                    */
    cosh_x = 0.5F*(exp_x + exp_minus_x);
    return cosh_x;
}
/*  End of rssringoccs_Float_Cosh.                                            */

/*  Double precision hyperbolic cosine (cosh equivalent).                     */
double rssringoccs_Double_Cosh(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double cosh_x, exp_x, exp_minus_x;

    /*  The definition of cosh(x) is [exp(x) + exp(-x)]/2, so return this. It *
     *  is computationally faster to compute exp(x) and then exp(-x) via the  *
     *  formula exp(-x) = 1/exp(x). This saves us from computing two          *
     *  exponentials at the cost of an extra division.                        */
    exp_x = rssringoccs_Double_Exp(x);
    exp_minus_x = 1.0/exp_x;

    /*  Compute cosh from the two exponentials and return.                    */
    cosh_x = 0.5*(exp_x + exp_minus_x);
    return cosh_x;
}
/*  End of rssringoccs_Double_Cosh.                                           */

/*  Long doubel precision hyperbolic cosine (coshl equivalent).               */
long double rssringoccs_LDouble_Cosh(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double cosh_x, exp_x, exp_minus_x;

    /*  The definition of cosh(x) is [exp(x) + exp(-x)]/2, so return this. It *
     *  is computationally faster to compute exp(x) and then exp(-x) via the  *
     *  formula exp(-x) = 1/exp(x). This saves us from computing two          *
     *  exponentials at the cost of an extra division.                        */
    exp_x = rssringoccs_LDouble_Exp(x);
    exp_minus_x = 1.0L/exp_x;

    /*  Compute cosh from the two exponentials and return.                    */
    cosh_x = 0.5L*(exp_x + exp_minus_x);
    return cosh_x;
}
/*  End of rssringoccs_LDouble_Cosh.                                          */
