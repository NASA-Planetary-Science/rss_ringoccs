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
 *      Contains the source code for the absolute value function.             *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Float_Abs:                                                *
 *      rssringoccs_Double_Abs:                                               *
 *      rssringoccs_LDouble_Abs:                                              *
 *  Purpose:                                                                  *
 *      Computes the absolute value of a real number.                         *
 *                   --                                                       *
 *                  |   x,  x > 0                                             *
 *          |x|  =  |  -x,  else                                              *
 *                   --                                                       *
 *  Arguments:                                                                *
 *      x (float/double/long double):                                         *
 *          A real number, the argument for |x|.                              *
 *  Output:                                                                   *
 *      abs_x (float/double/long double):                                     *
 *          The absolute value of x.                                          *
 *  Method:                                                                   *
 *      This uses a simple if-then statement to check if the input is         *
 *      positive or not, returning x for non-negative and -x otherwise.       *
 *  Notes:                                                                    *
 *      fabsf and fabsl are not provided in C89/C90 implementations of the    *
 *      language, and instead type conversions are made in the fabs function. *
 *      Since the absolute value function is very simple, we simply provide   *
 *      the algorithms here rather than pass the arguments to fabs, fabsf, or *
 *      fabsfl. There is essentially no time difference. Using gcc with -O2   *
 *      optimization on an array of 1 million random elements in [-1, 1] gave *
 *      the following times (in seconds):                                     *
 *          fabs:         0.003328                                            *
 *          rss_ringoccs: 0.003743                                            *
 *      -O3 optimization gives:                                               *
 *          fabs:         0.003409                                            *
 *          rss_ringoccs: 0.003493                                            *
 *      So, no real difference. These times were computed on an iMac 2017     *
 *      3.4GHz quad-core running MacOS Catalina 10.15.7. Converting a long    *
 *      double to a double may lose precision, hence the reason we provide    *
 *      this simple code.                                                     *
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

/*  Single precision absolute value function (fabsf equivalent).              */
float rssringoccs_Float_Abs(float x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    float abs_x;

    /*  If x is positive return it, otherwise return its negative.            */
    if (x >= 0.0F)
        abs_x = x;
    else
        abs_x = -x;

    return abs_x;
}
/*  End of rssringoccs_Float_Abs.                                             */

/*  Double precision absolute value function (fabs equivalent).               */
double rssringoccs_Double_Abs(double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double abs_x;

    /*  If x is positive return it, otherwise return its negative.            */
    if (x >= 0.0)
        abs_x = x;
    else
        abs_x = -x;

    return abs_x;
}
/*  End of rssringoccs_Double_Abs.                                            */

/*  Long double precision absolute value function (fabsl equivalent).         */
long double rssringoccs_LDouble_Abs(long double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    long double abs_x;

    /*  If x is positive return it, otherwise return its negative.            */
    if (x >= 0.0L)
        abs_x = x;
    else
        abs_x = -x;

    return abs_x;
}
/*  End of rssringoccs_LDouble_Abs.                                           */
