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

/*  Library for timing computations.                                          */
#include <time.h>

/*  Needed for printing the outputs.                                          */
#include <stdio.h>

/*  Needed for malloc and exit.                                               */
#include <stdlib.h>

/*  Needed for fabs, fabsf, and fabsl.                                        */
#include <math.h>

/*  Prototypes for these functions found here.                                */
#include "rss_ringoccs_compare_funcs.h"

/*  rss_ringoccs complex variables defined here.                              */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  C99 standard complex library.                                             */
#include <complex.h>

/*  Routine for comparing complex valued functions of a double precision      *
 *  complex variable.                                                         */
void
rssringoccs_Compare_CDouble_Funcs(
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

    /*  Declare variables for computing |z0-z1|.                              */
    double x_c, y_c, x_s, y_s, dx, dy;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    double max_err = 0.0;
    double temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    rssringoccs_ComplexDouble **z0;
    double _Complex **z1;

    /*  Declare a dummy variable for indexing.                                */
    unsigned long m, n;

    /*  Declare variables for computing computation time.                     */
    clock_t t1, t2;

    /*  Error check to make sure the user provided valid inputs.              */
    if (start >= end)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_CDouble_Funcs\n\n"
             "start is greater than or equal to end.\n"
             "Abort computation.\n");
        exit(0);
    }
    else if (N == 0)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_CDouble_Funcs\n\n"
             "Input sample size is zero. Aborting computation.\n");
        exit(0);
    }

    /*  We'll increment evenly throughout the region.                         */
    ds = (end - start) / N;

    /*  Allocate memory for the two pointers we've declared.                  */
    z0 = malloc(sizeof(*z0) * N);
    z1 = malloc(sizeof(*z1) * N);

    /*  And set each entry of z0 and z1 to point to blocks of memory for N    *
     *  complex doubles.                                                      */
    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }
    /*  Be sure to free when done!                                            */

    /*  Set y to the starting value and grab the current time.                */
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        /*  Reset x to the starting value.                                    */
        x = start;
        for (n=0; n<N; ++n)
        {
            /*  Story the value of the (n, m) entry in z0[n][m].              */
            z0[n][m] = f0(rssringoccs_CDouble_Rect(x, y));

            /*  Increment x.                                                  */
            x += ds;
        }

        /*  Increment y. x will be reset above for the next iteration.        */
        y += ds;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (m=0; m<N; ++m)
    {
        /*  Reset x to the start.                                             */
        x = start;
        for (n=0; n<N; ++n)
        {
            /*  Store the value of the (n, m) entry in z1[n][m].              */
            z1[n][m] = f1(x + (double _Complex)_Complex_I * y);

            /*  Increment x.                                                  */
            x += ds;
        }

        /*  Increment y. x will be reset above for the next iteration.        */
        y += ds;
    }

    /*  Grab the time again.                                                  */
    t2 = clock();

    /*  Print out how long it took for f1 to compute.                         */
    printf("%s: %f\n", f1_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  NOTE:                                                                 *
     *      Without the following comparison of the two pointers y0 and y1,   *
     *      some compilers may see the above computations as redundant with   *
     *      optimization on, and skip them. The resulting times will be close *
     *      to zero for both f0 and f1.                                       */

    /*  Compute the maximum absolute error between f0 and f1.                 */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            /*  Extract the real and imaginary part from z0.                  */
            x_s = rssringoccs_CDouble_Real_Part(z0[n][m]);
            y_s = rssringoccs_CDouble_Imag_Part(z0[n][m]);

            /*  Extract the real and imaginary part from z1.                  */
            x_c = creal(z1[n][m]);
            y_c = cimag(z1[n][m]);

            /*  We'll use the Pythagorean formula to compute |z0-z1|. Compute *
             *  the difference in the real and imaginary parts of z0 and z1.  */
            dx = x_s - x_c;
            dy = y_s - y_c;

            /*  We'll use the standard library function to check the error.   */
            temp = sqrt(dx*dx + dy*dy);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |z0-z1|.                                    */

    /*  Print out the error to 16 decimals (assumes 64-bit precision).        */
    printf("Max Error: %.16f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    /*  And free the pointers to the pointers.                                */
    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_CDouble_Funcs.                                 */
