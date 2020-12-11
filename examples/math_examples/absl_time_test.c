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
 *      Provide tests for the accuracy and efficiency of rss_ringoccs         *
 *      absolute value function at long double precision compared to the one  *
 *      provided by the C99 standard in math.h.                               *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  The absolute value functions are found here.                              */
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <math.h>

/*  Library for timing computations.                                          */
#include <time.h>

/*  Needed for printing the outputs.                                          */
#include <stdio.h>

/*  Needed for malloc.                                                        */
#include <stdlib.h>

int main(void)
{
    /*  Declare necessary variables.                                          */
    long double max_err, temp, start, end, x, dx;
    long double *y0, *y1;
    unsigned int n, N, ind;
    clock_t t1, t2;

    /*  We'll do our time test with 100 million points.                       */
    N = 1e8;

    /*  We'll have the variable range from -100 to 100.                       */
    start = 100.0L;
    end = 100.0L;

    /*  Set the initial values for ind and max_err to zero.                   */
    ind = 0;
    max_err = 0.0L;

    /*  And we'll increment evenly throughout the region.                     */
    dx = (end - start) / N;

    /*  Allocate memory for the two pointers we've declared.                  */
    y0 = malloc(sizeof(*y0) * N);
    y1 = malloc(sizeof(*y1) * N);

    /*  Set x to the starting value and grab the current time.                */
    x = -start;
    t1 = clock();

    /*  Perform the calculation for the C99 standard library function fabsl.  */
    for (n=0; n<N; ++n)
    {
        y0[n] = fabsl(x);
        x += dx;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to second, use the macro CLOCKS_PER_SEC    *
     *  provided in time.h.                                                   */
    printf("C99:          %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the rss_ringoccs function.                */
    x = -start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using rssringoccs_LDouble_Abs instead.        */
    for (n=0; n<N; ++n)
    {
        y1[n] = rssringoccs_LDouble_Abs(x);
        x += dx;
    }

    /*  Grab the time again.                                                  */
    t2 = clock();

    /*  Print out how long it took for rss_ringoccs to compute.               */
    printf("rss_ringoccs: %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Compute the maximum absolute error between rss_ringoccs and C99.      */
    for (n=0; n<N; ++n)
    {
        /*  We'll use the standard library function to check the error.       */
        temp = fabsl(y0[n] - y1[n]);

        /*  Check if the error got larger and set max_err accordingly.        */
        if (max_err < temp)
        {
            max_err = temp;
            ind = n;
        }
    }
    /*  End of for loop computing |y0-y1|.                                    */

    /*  Print out the error to 24 decimals (assumes 96-bit precision).        */
    printf("Max Error: %.24Lf\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    free(y0);
    free(y1);
    return 0;
}

/*  This was compiled on an iMac 2017 3.4GHz quad-core running MacOS Catalina *
 *  10.15.7. It produces the following times:                                 *
 *      NOTE: On MacOS gcc is aliased to LLVM's clang:                        *
 *          gcc --version                                                     *
 *          Apple clang version 12.0.0 (clang-1200.0.32.27)                   *
 *      This is NOT the regular gcc from GNU. To use this on apple devices    *
 *      requires homebrew.                                                    *
 *                                                                            *
 *  -O3 optimization:                                                         *
 *      gcc -Wconversion -O3 -Wall -Wextra -Wpedantic                         *
 *          -pedantic absl_time_test.c -o test -lrssringoccs                  *
 *      C99:          0.667850                                                *
 *      rss_ringoccs: 1.083205                                                *
 *      Max Error: 0.000000000000000000000000                                 *
 *  No optimization:                                                          *
 *      gcc absl_time_test.c -o test -lrssringoccs                            *
 *      C99:          1.154483                                                *
 *      rss_ringoccs: 1.566479                                                *
 *      Max Error: 0.000000000000000000000000                                 */
