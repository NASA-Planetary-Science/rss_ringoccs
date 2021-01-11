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
 *  Date:       December 29, 2020                                             *
 ******************************************************************************/

/*  Library for timing computations.                                          */
#include <time.h>

/*  Needed for printing the outputs.                                          */
#include <stdio.h>

/*  Needed for malloc.                                                        */
#include <stdlib.h>

/*  Needed for fabs, fabsf, and fabsl.                                        */
#include <math.h>

/*  Prototypes for these functions found here.                                */
#include "rss_ringoccs_compare_funcs.h"

/*  Routine for comparing two real valued functions at long double precision. */
void
rssringoccs_Compare_LDouble_Funcs(
    const char *f0_name,
    long double (*f0)(long double),
    const char *f1_name,
    long double (*f1)(long double),
    const long double start,
    const long double end,
    const unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    long double x, dx;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    long double max_err = 0.0L;
    long double temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    long double *y0, *y1;

    /*  Declare a dummy variable for indexing.                                */
    unsigned long n;

    /*  Declare variables for computing computation time.                     */
    clock_t t1, t2;

    /*  Error check to make sure the user provided valid inputs.              */
    if (start >= end)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_LDouble_Funcs\n\n"
             "start is greater than or equal to end.\n"
             "Abort computation.\n");
        return;
    }
    else if (N == 0)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_LDouble_Funcs\n\n"
             "Input sample size is zero. Aborting computation.\n");
        return;
    }

    /*  We'll increment evenly throughout the region.                         */
    dx = (end - start) / N;

    /*  Allocate memory for the two pointers we've declared.                  */
    y0 = malloc(sizeof(*y0) * N);
    y1 = malloc(sizeof(*y1) * N);

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (n=0; n<N; ++n)
    {
        y0[n] = f0(x);
        x += dx;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    x = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (n=0; n<N; ++n)
    {
        y1[n] = f1(x);
        x += dx;
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
    for (n=0; n<N; ++n)
    {
        /*  We'll use the standard library function to check the error.       */
        temp = fabsl(y0[n] - y1[n]);

        /*  Check if the error got larger and set max_err accordingly.        */
        if (max_err < temp)
            max_err = temp;
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 24 decimals (assumes 96-bit precision).        */
    printf("Max Error: %.24Lf\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    free(y0);
    free(y1);
}
/*  End of rssringoccs_Compare_LDouble_Funcs.                                 */
