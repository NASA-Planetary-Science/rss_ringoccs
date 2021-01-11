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
 *                 rss_ringoccs_sorted_linear_interpolation                   *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for linear interpolation of sorted data.     *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Float_Sorted_Interp1d:                                    *
 *      rssringoccs_Double_Sorted_Interp1d:                                   *
 *      rssringoccs_LDouble_Sorted_Interp1d:                                  *
 *  Purpose:                                                                  *
 *      Computes the absolute value, or modulus, of a complex number:         *
 *                                                                            *
 *          |z| = |x + iy| = sqrt(x^2 + y^2)                                  *
 *                                                                            *
 *  Arguments:                                                                *
 *      x (float *, double *, long double *):                                 *
 *          A sorted array of real numbers that are strictly                  *
 *          monotonically increasing.                                         *
 *      y (float *, double *, long double *):                                 *
 *          The data points corresponding to x.                               *
 *      N (unsigned long):                                                    *
 *          The number of elements of x and y.                                *
 *      x_new (float *, double *, long double *):                             *
 *          The new data points. Must lie between min(x) and max(x).          *
 *      y_new (float *, double *, long double *):                             *
 *          The interpolated data corresponding to x_new.                     *
 *      N_new (unsigned long):                                                *
 *          The number of elements of x_new and y_new.                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Method:                                                                   *
 *      Find the least value x[n] such x_new[m] < x[n] and then perform a     *
 *      linear interpolation with x[n-1] and x[n] via the formula:            *
 *                                                                            *
 *                               y[n] - y[n-1]                                *
 *          y_new[m] = y[n-1] + --------------- * (x_new[m] - x[n-1])         *
 *                               x[n] - x[n-1]                                *
 *                                                                            *
 *  NOTES:                                                                    *
 *      No error checks are made on whether or not the pointers are NULL or   *
 *      if there are N and N_new elements to x, y, and x_new, y_new,          *
 *      respectively. Improper use of this routine may cause a segmentation   *
 *      fault.                                                                *
 *                                                                            *
 *      It is assumed the data is ordered and strictly monotonically          *
 *      increasing. If these assumptions are not true the outcome will be     *
 *      meaningless and some values may return NaN.                           *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_interpolate.h:                                           *
 *          Header where the function prototypes are defined.                 *
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
 *  Date:       December 30, 2020                                             *
 ******************************************************************************/

/*  rssringoccs_NaN is defined here.                                          */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  And the function prototypes are found here.                               */
#include <rss_ringoccs/include/rss_ringoccs_interpolate.h>

/*  Single precision linear interpolation of sorted data.                     */
void
rssringoccs_Float_Sorted_Interp1d(float *x,
                                  float *y,
                                  unsigned long N,
                                  float *x_new,
                                  float *y_new,
                                  unsigned long N_new)
{
    /*  Declare two variables for indexing the interpolated and raw data.     */
    unsigned long m, n;

    /*  And declare a variable for computing the slope for the interpolation. */
    float slope;

    /*  Loop over the entries of the interpolated pointers and compute.       */
    for (m=0; m<N_new; ++m)
    {
        /*  Reset n to zero.                                                  */
        n = 0;

        /*  If x_new[m] falls outside of the bounds of x, return NaN.         */
        if ((x_new[m] < x[0]) || (x_new[m] > x[N-1]))
            y_new[m] = rssringoccs_NaN_F;

        /*  Handle the case of x_new[m] = x[N-1].                             */
        else if (x_new[m] == x[N-1])
            y_new[m] = y[N-1];

        /*  Handle the special case of x_new[m] = x[0].                       */
        else if (x_new[m] == x[0])
            y_new[m] = y[0];

        /*  And finally, handle the general case.                             */
        else
        {
            /*  Find the smallest index n such that x[n] > x_new[m].          */
            while (x[n] <= x_new[m])
                n++;

            /*  Use this index to compute the linear interpolation.           */
            slope = (y[n] - y[n-1]) / (x[n] - x[n-1]);
            y_new[m] = y[n-1] + slope * (x_new[m] - x[n-1]);
        }
    }
    /*  End of for loop computing y_new[m].                                   */
}
/*  End of rssringoccs_Float_Sorted_Interp1d.                                 */

/*  Double precision linear interpolation of sorted data.                     */
void
rssringoccs_Double_Sorted_Interp1d(double *x,
                                   double *y,
                                   unsigned long N,
                                   double *x_new,
                                   double *y_new,
                                   unsigned long N_new)
{
    /*  Declare two variables for indexing the interpolated and raw data.     */
    unsigned long m, n;

    /*  And declare a variable for computing the slope for the interpolation. */
    double slope;

    /*  Loop over the entries of the interpolated pointers and compute.       */
    for (m=0; m<N_new; ++m)
    {
        /*  Reset n to zero.                                                  */
        n = 0;

        /*  If x_new[m] falls outside of the bounds of x, return NaN.         */
        if ((x_new[m] < x[0]) || (x_new[m] > x[N-1]))
            y_new[m] = rssringoccs_NaN;

        /*  Handle the case of x_new[m] = x[N-1].                             */
        else if (x_new[m] == x[N-1])
            y_new[m] = y[N-1];

        /*  Handle the special case of x_new[m] = x[0].                       */
        else if (x_new[m] == x[0])
            y_new[m] = y[0];

        /*  And finally, handle the general case.                             */
        else
        {
            /*  Find the smallest index n such that x[n] > x_new[m].          */
            while (x[n] <= x_new[m])
                n++;

            /*  Use this index to compute the linear interpolation.           */
            slope = (y[n] - y[n-1]) / (x[n] - x[n-1]);
            y_new[m] = y[n-1] + slope * (x_new[m] - x[n-1]);
        }
    }
    /*  End of for loop computing y_new[m].                                   */
}
/*  End of rssringoccs_Double_Sorted_Interp1d.                                */

/*  Long double precision linear interpolation of sorted data.                */
void
rssringoccs_LDouble_Sorted_Interp1d(long double *x,
                                    long double *y,
                                    unsigned long N,
                                    long double *x_new,
                                    long double *y_new,
                                    unsigned long N_new)
{
    /*  Declare two variables for indexing the interpolated and raw data.     */
    unsigned long m, n;

    /*  And declare a variable for computing the slope for the interpolation. */
    long double slope;

    /*  Loop over the entries of the interpolated pointers and compute.       */
    for (m=0; m<N_new; ++m)
    {
        /*  Reset n to zero.                                                  */
        n = 0;

        /*  If x_new[m] falls outside of the bounds of x, return NaN.         */
        if ((x_new[m] < x[0]) || (x_new[m] > x[N-1]))
            y_new[m] = rssringoccs_NaN_L;

        /*  Handle the case of x_new[m] = x[N-1].                             */
        else if (x_new[m] == x[N-1])
            y_new[m] = y[N-1];

        /*  Handle the special case of x_new[m] = x[0].                       */
        else if (x_new[m] == x[0])
            y_new[m] = y[0];

        /*  And finally, handle the general case.                             */
        else
        {
            /*  Find the smallest index n such that x[n] > x_new[m].          */
            while (x[n] <= x_new[m])
                n++;

            /*  Use this index to compute the linear interpolation.           */
            slope = (y[n] - y[n-1]) / (x[n] - x[n-1]);
            y_new[m] = y[n-1] + slope * (x_new[m] - x[n-1]);
        }
    }
    /*  End of for loop computing y_new[m].                                   */
}
/*  End of rssringoccs_LDouble_Sorted_Interp1d.                               */
