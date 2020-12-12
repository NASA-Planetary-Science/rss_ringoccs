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
 *      absolute value function at single precision compared to the one       *
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

/*  Routine for comparing fabsf with rssringoccs_Float_Abs.                   */
int main(void)
{
    /*  Set the start and end for the values we're testing.                   */
    float start = -100.0;
    float end   =  100.0;

    /*  Declare variables for sampling the region [start, end].               */
    float x, dx;

    /*  Declare variables for computing the maximum difference between fabsf  *
     *  and rssringoccs_Float_Abs.                                            */
    float max_err = 0.0;
    float temp;

    /*  Declare two pointers to represent arrays for fabsf(x) and             *
     *  rssringoccs_Float_Abs(x), respectively.                               */
    float *y0, *y1;

    /*  Declare a dummy variable for indexing and a variable for the number   *
     *  of points we're sampling in the range [start, end].                   */
    unsigned int n;
    unsigned int N = 1e8;

    /*  Declare variables for computing computation time.                     */
    clock_t t1, t2;

    /*  We'll increment evenly throughout the region.                         */
    dx = (end - start) / N;

    /*  Allocate memory for the two pointers we've declared.                  */
    y0 = malloc(sizeof(*y0) * N);
    y1 = malloc(sizeof(*y1) * N);

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    t1 = clock();

    /*  Perform the calculation for the C99 standard library function fabsf.  */
    for (n=0; n<N; ++n)
    {
        y0[n] = fabsf(x);
        x += dx;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("C99:          %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the rss_ringoccs function.                */
    x = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using rssringoccs_Float_Abs instead of fabsf. */
    for (n=0; n<N; ++n)
    {
        y1[n] = rssringoccs_Float_Abs(x);
        x += dx;
    }

    /*  Grab the time again.                                                  */
    t2 = clock();

    /*  Print out how long it took for rss_ringoccs to compute.               */
    printf("rss_ringoccs: %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  NOTE:                                                                 *
     *      Without the following comparison of the two pointers y0 and y1,   *
     *      some compilers may see the above computations as redundant with   *
     *      optimization on, and skip them. The resulting times will be close *
     *      to zero for both fabsf and rssringoccs_Float_Abs.                 */

    /*  Compute the maximum absolute error between rss_ringoccs and C99.      */
    for (n=0; n<N; ++n)
    {
        /*  We'll use the standard library function to check the error.       */
        temp = fabsf(y0[n] - y1[n]);

        /*  Check if the error got larger and set max_err accordingly.        */
        if (max_err < temp)
            max_err = temp;
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.8f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    free(y0);
    free(y1);
    return 0;
}
/*  End of main.                                                              */

/*  This was compiled on an iMac 2017 3.4GHz quad-core running MacOS Catalina *
 *  10.15.7. It produces the following times:                                 *
 *      NOTE: On MacOS gcc is aliased to LLVM's clang:                        *
 *          gcc --version                                                     *
 *          Apple clang version 12.0.0 (clang-1200.0.32.27)                   *
 *      This is NOT the regular gcc from GNU. To use this on apple devices    *
 *      requires homebrew.                                                    *
 *                                                                            *
 *  No optimization:                                                          *
 *      gcc -Wall -Wextra -Wconversion -Wpedantic -pedantic -pedantic-errors  *
 *          -std=c89 -ansi absf_time_test.c -o test -lrssringoccs             *
 *      C99:          0.374127                                                *
 *      rss_ringoccs: 0.405311                                                *
 *      Max Error: 0.00000000                                                 *
 *  -O3 optimization:                                                         *
 *      gcc -Wall -Wextra -Wconversion -Wpedantic -pedantic -pedantic-errors  *
 *          -O3 -std=c89 -ansi absf_time_test.c -o test -lrssringoccs         *
 *      C99:          0.129045                                                *
 *      rss_ringoccs: 0.290840                                                *
 *      Max Error: 0.00000000                                                 */
