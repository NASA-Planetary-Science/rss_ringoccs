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
 *      absolute value function compared to the one provided by the C99       *
 *      standard in math.h.                                                   *
 *  NOTE:                                                                     *
 *      If rss_ringoccs was built with C99 math.h support, i.e. the macro     *
 *      __RSS_RINGOCCS_USING_C99_MATH_H__ was set to 1 in                     *
 *      rss_ringoccs_config.h, then this test is redundant since fabs and     *
 *      rssringoccs_Double_Abs are the same thing.                            *
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
    double max_err, temp, start, end, x, dx;
    double *y0, *y1;
    unsigned int n, N, ind;
    clock_t t1, t2;

    /*  We'll do our time test with 100 million points.                       */
    N = 1e8;

    /*  We'll have the variable range from -100 to 100.                       */
    start = 100.0;
    end = 100.0;

    /*  And we'll increment evenly throughout the region.                     */
    dx = (end - start) / N;

    /*  Set the initial values for ind and max_err to zero.                   */
    ind = 0;
    max_err = 0.0;

    /*  Allocate memory for the two pointers we've declared.                  */
    y0 = malloc(sizeof(*y0) * N);
    y1 = malloc(sizeof(*y1) * N);

    /*  Set x to the starting value and grab the current time.                */
    x = -start;
    t1 = clock();

    /*  Perform the calculation for the C99 standard library function fabs.   */
    for (n=0; n<N; ++n)
    {
        y0[n] = fabs(x);
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

    /*  Perform the computation using rssringoccs_Double_Abs instead of fabs. */
    for (n=0; n<N; ++n)
    {
        y1[n] = rssringoccs_Double_Abs(x);
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
     *      to zero for both fabs and rssringoccs_Double_Abs.                 */

    /*  Compute the maximum absolute error between rss_ringoccs and C99.      */
    for (n=0; n<N; ++n)
    {
        /*  We'll use the standard library function to check the error.       */
        temp = fabs(y0[n] - y1[n]);

        /*  Check if the error got larger and set max_err accordingly.        */
        if (max_err < temp)
        {
            max_err = temp;
            ind = n;
        }
    }
    /*  End of for loop computing |y0-y1|.                                    */

    /*  Print out the error to 16 decimals (assumes 64-bit precision).        */
    printf("Max Error: %.16f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    free(y0);
    free(y1);
    return 0;
}
/*  End of main.                                                              */

/*  This was compiled with various options on an iMac 2017 3.4GHz quad-core   *
 *  running MacOS Catalina 10.15.7. It produces the following times:          *
 *      NOTE: On MacOS gcc is aliased to LLVM's clang:                        *
 *          gcc --version                                                     *
 *          Apple clang version 12.0.0 (clang-1200.0.32.27)                   *
 *      This is NOT the regular gcc from GNU. To use gcc on apple devices     *
 *      requires homebrew.                                                    *
 *  c89 option, no optimization:                                              *
 *      gcc -Wall -Wextra -Wpedantic -pedantic -std=c89                       *
 *              -ansi abs_time_test.c -o test -lrssringoccs                   *
 *      C99:          0.513468                                                *
 *      rss_ringoccs: 0.537242                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  c89 option, O2 optimization.                                              *
 *      gcc -O2 -Wall -Wextra -Wpedantic -pedantic -std=c89                   *
 *              -ansi abs_time_test.c -o test -lrssringoccs                   *
 *      C99:          0.264355                                                *
 *      rss_ringoccs: 0.449514                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  c89 option, O3 optimization:                                              *
 *      gcc -O2 -Wall -Wextra -Wpedantic -pedantic -std=c89                   *
 *              -ansi abs_time_test.c -o test -lrssringoccs                   *
 *      C99:          0.268884                                                *
 *      rss_ringoccs: 0.446283                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  c99 option, no optimization.                                              *
 *       gcc  -Wall -Wextra -Wpedantic -pedantic                              *
 *              -std=c99 abs_time_test.c -o test -lrssringoccs                *
 *      C99:          0.500101                                                *
 *      rss_ringoccs: 0.528889                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  c99 option, -O2 optimization:                                             *
 *      gcc -O2 -Wall -Wextra -Wpedantic -pedantic                            *
 *              -std=c99 abs_time_test.c -o test -lrssringoccs                *
 *      C99:          0.264327                                                *
 *      rss_ringoccs: 0.446251                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  c99 option, -O3 optimization:                                             *
 *      gcc -O3 -Wall -Wextra -Wpedantic -pedantic                            *
 *              -std=c99 abs_time_test.c -o test -lrssringoccs                *
 *      C99:          0.268204                                                *
 *      rss_ringoccs: 0.448208                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  Default, -O3 optimization.                                                *
 *      gcc -O3 abs_time_test.c -o test -lrssringoccs                         *
 *      C99:          0.284114                                                *
 *      rss_ringoccs: 0.465479                                                *
 *      Max Error: 0.0000000000000000                                         *
 *  NOTE:                                                                     *
 *      These times will differ on different devices and on different         *
 *      compilers.                                                            */
