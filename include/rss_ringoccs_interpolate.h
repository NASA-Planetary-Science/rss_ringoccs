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
 *                         rss_ringoccs_interpolate                           *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Provide routines for interpolating data.                              *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 30, 2020                                             *
 ******************************************************************************/

/*  Include guard to prevent including this file twice.                       */
#ifndef __RSS_RINGOCCS_INTERPOLATE_H__
#define __RSS_RINGOCCS_INTERPOLATE_H__

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Double_Sorted_Interp1d                                    *
 *  Purpose:                                                                  *
 *      Interpolate values of (x, y) to (x_new, y_new).                       *
 *  Arguments:                                                                *
 *      x (double *):                                                         *
 *          A sorted array of real numbers that are strictly                  *
 *          monotonically increasing.                                         *
 *      y (double *):                                                         *
 *          The data points corresponding to x.                               *
 *      N (unsigned long):                                                    *
 *          The number of elements of x and y.                                *
 *      x_new (double *):                                                     *
 *          The new data points. Must lie between min(x) and max(x).          *
 *      y_new (double *):                                                     *
 *          The interpolated data corresponding to x_new.                     *
 *      N_new (unsigned long):                                                *
 *          The number of elements of x_new and y_new.                        *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  NOTES:                                                                    *
 *      It is assumed the pointer x has it's data sorted and strictly         *
 *      monotonically increasing. That is, x[n] < x[n+1] for all valid n.     *
 ******************************************************************************/
extern void
rssringoccs_Float_Sorted_Interp1d(float *x,
                                  float *y,
                                  unsigned long N,
                                  float *x_new,
                                  float *y_new,
                                  unsigned long N_new);

extern void
rssringoccs_Double_Sorted_Interp1d(double *x,
                                   double *y,
                                   unsigned long N,
                                   double *x_new,
                                   double *y_new,
                                   unsigned long N_new);

extern void
rssringoccs_LDouble_Sorted_Interp1d(long double *x,
                                    long double *y,
                                    unsigned long N,
                                    long double *x_new,
                                    long double *y_new,
                                    unsigned long N_new);

#endif
/*  End of include guard.                                                     */
