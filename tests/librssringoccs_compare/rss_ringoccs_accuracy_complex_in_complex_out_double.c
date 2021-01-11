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

/*  Routine for comparing two complex valued functions at single precision.   */
void
rssringoccs_Accuracy_CDouble_Funcs(
    const char *f0_name,
    rssringoccs_ComplexDouble (*f0)(rssringoccs_ComplexDouble),
    const char *f1_name,
    double _Complex (*f1)(double _Complex),
    const double start,
    const double end,
    const unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    double x, y, ds;

    /*  Declare variables for our complex values.                             */
    rssringoccs_ComplexDouble z0;
    double _Complex z1;

    /*  Declare variables for computing |z0-z1|/|z1|.                         */
    double x_c, y_c, x_s, y_s, dx, dy, x_w, y_w, f0_x, f0_y, f1_x, f1_y;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    double max_err = 0.0;
    double temp;

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

    /*  Initialize the variables for the "worst" error in the plane.          */
    x_w = start;
    y_w = start;
    z0 = f0(rssringoccs_CDouble_Rect(x_w, y_w));
    z1 = f1(x_w + (double _Complex)_Complex_I*y_w);
    f0_x = rssringoccs_CDouble_Real_Part(z0);
    f0_y = rssringoccs_CDouble_Imag_Part(z0);
    f1_x = creal(z1);
    f1_y = cimag(z1);

    /*  Compute the maximum absolute error between f0 and f1.                 */
    for (m=0; m<N; ++m)
    {
        /*  Reset x to its starting value.                                    */
        x = start;
        for (n=0; n<N; ++n)
        {
            z0 = f0(rssringoccs_CDouble_Rect(x, y));
            z1 = f1(x + (double _Complex)_Complex_I*y);

            /*  Extract the real and imaginary part from z0.                  */
            x_s = rssringoccs_CDouble_Real_Part(z0);
            y_s = rssringoccs_CDouble_Imag_Part(z0);

            /*  Extract the real and imaginary part from z1.                  */
            x_c = creal(z1);
            y_c = cimag(z1);

            /*  Compute the difference of the real and imaginary parts of     *
             *  z0 and z1. We'll use this together with the Pythagorean       *
             *  formula to compute |z0-z1|.                                   */
            dx = x_s - x_c;
            dy = y_s - y_c;

            /*  We'll use the standard library function to check the error.   */
            temp = sqrt((dx*dx + dy*dy) / (x_c*x_c + y_c*y_c));

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
            {
                max_err = temp;
                x_w = x;
                y_w = y;
                f0_x = x_s;
                f0_y = y_s;
                f1_x = x_c;
                f1_y = y_c;
            }

            /*  Increment x.                                                  */
            x += ds;
        }

        /*  Increment y.                                                      */
        y += ds;
    }
    /*  End of for-loop computing |z0-z1|.                                    */

    /*  Print out the error to 16 decimals (assumes 64-bit precision).        */
    printf("\n%s vs %s:\n", f0_name, f1_name);
    printf("\tRegion of Complex Plane Sampled:\n\t\t [%f, %f]x[%f, %f]\n",
           start, end, start, end);
    printf("\tNumber of Points Sampled:\n\t\t%lux%lu\n", N, N);
    printf("\tMaximum Relative Error: %.16f\n", max_err);
    printf("\tPoint of Highest Error:\n\t\t(%.16f, %.16f)\n", x_w, y_w);
    printf("\t%s Value:\n\t\t%.16f + %.16fi\n", f0_name, f0_x, f0_y);
    printf("\t%s Value:\n\t\t%.16f + %.16fi\n", f1_name, f1_x, f1_y);
}
/*  End of rssringoccs_Accuracy_CDouble_Funcs.                                */
