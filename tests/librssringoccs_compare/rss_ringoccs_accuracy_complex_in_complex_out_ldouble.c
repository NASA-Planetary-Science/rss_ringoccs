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
 *  Purpose:                                                                  *
 *      Provide functions for comparing the accuracy and efficiency of        *
 *      functions in rss_ringoccs as opposed to other libraries.              *
 *  NOTE:                                                                     *
 *      librssringoccs does not have any dependencies and will compile on any *
 *      compiler capable of handling C89/C90 or C99 compliant code. The tests *
 *      using these functions use external libraries to compare the results   *
 *      of rss_ringoccs with others. To run these tests requires having these *
 *      libraries available. These tests are NOT required to use rss_ringoccs *
 *      and are mainly for internal use.                                      *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 12, 2020                                             *
 ******************************************************************************/

/*  Needed for printing the outputs.                                          */
#include <stdio.h>

/*  Needed for fabs, fabsf, and fabsl.                                        */
#include <math.h>

/*  Prototypes for these functions found here.                                */
#include "rss_ringoccs_compare_funcs.h"
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <complex.h>

extern void
rssringoccs_Accuracy_CLDouble_Funcs(
    const char *f0_name,
    rssringoccs_ComplexLongDouble (*f0)(rssringoccs_ComplexLongDouble),
    const char *f1_name,
    long double _Complex (*f1)(long double _Complex),
    const long double start,
    const long double end,
    const unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    long double x, y, ds;

    /*  Declare variables for our complex values.                             */
    rssringoccs_ComplexLongDouble z0;
    long double _Complex z1;

    /*  Declare variables for computing |z0-z1|/|z1|.                         */
    long double x_c, y_c, x_s, y_s, dx, dy;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    long double max_err = 0.0L;
    long double temp;

    /*  Declare a dummy variable for indexing.                                */
    unsigned long m, n;

    /*  Error check to make sure the user provided valid inputs.              */
    if (start >= end)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Accuracy_CDouble_Funcs\n\n"
             "start is greater than or equal to end.\n"
             "Abort computation.\n");
        return;
    }
    else if (N == 0)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Accuracy_CDouble_Funcs\n\n"
             "Input sample size is zero. Aborting computation.\n");
        return;
    }

    /*  We'll increment evenly throughout the region.                         */
    ds = (end - start) / N;

    /*  Set y to it's initial value.                                          */
    y = start;

    /*  Compute the maximum absolute error between f0 and f1.                 */
    for (m=0; m<N; ++m)
    {
        /*  Reset x to its starting value.                                    */
        x = start;
        for (n=0; n<N; ++n)
        {
            z0 = f0(rssringoccs_CLDouble_Rect(x, y));
            z1 = f1(x + (long double _Complex)_Complex_I*y);

            /*  Extract the real and imaginary part from z0.                  */
            x_s = rssringoccs_CLDouble_Real_Part(z0);
            y_s = rssringoccs_CLDouble_Imag_Part(z0);

            /*  Extract the real and imaginary part from z1.                  */
            x_c = creall(z1);
            y_c = cimagl(z1);

            /*  Compute the difference of the real and imaginary parts of     *
             *  z0 and z1. We'll use this together with the Pythagorean       *
             *  formula to compute |z0-z1|.                                   */
            dx = x_s - x_c;
            dy = y_s - y_c;

            /*  We'll use the standard library function to check the error.   */
            temp = sqrtl((dx*dx + dy*dy) / (x_c*x_c + y_c*y_c));

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;

            /*  Increment x.                                                  */
            x += ds;
        }

        /*  Increment y.                                                      */
        y += ds;
    }
    /*  End of for-loop computing |z0-z1|.                                    */

    /*  Print out the error to 24 decimals (assumes 96-bit precision).        */
    printf("%s vs %s:\n", f0_name, f1_name);
    printf("\tRegion of Complex Plane Sampled: [%Lf, %Lf]x[%Lf, %Lf]\n",
           start, end, start, end);
    printf("\tNumber of Points Sampled: %lu\n", N);
    printf("\tMaximum Relative Error: %.24Lf\n", max_err);
}
/*  End of rssringoccs_Accuracy_CDouble_Funcs.                                */
