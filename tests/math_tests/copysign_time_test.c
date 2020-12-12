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
 *      copysign function compared to the C99 version.                        *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  The copysign functions are found here.                                    */
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
    double max_err, temp, start, end, x, dx, y, dy;
    double **z0, **z1;
    unsigned int m, n, N;
    clock_t t1, t2;

    /*  We'll do our time test with a square of 1000x1000 points.             */
    N = 1e4;

    /*  We'll have the variables range from -10 to 10.                        */
    start = -10.0;
    end = 10.0;

    /*  And we'll increment evenly throughout the region.                     */
    dx = (end - start) / N;
    dy = dx;

    /*  Set the initial value for max_err to zero.                            */
    max_err = 0.0;

    /*  Allocate memory for the two pointers we've declared.                  */
    z0 = malloc(sizeof(*z0) * N);
    z1 = malloc(sizeof(*z1) * N);

    /*  Allocate memory for each entry of the z0 and z1 pointers.             */
    for (n=0; n<N; ++n)
    {
        z0[n] = malloc(sizeof(*z0[n]) * N);
        z1[n] = malloc(sizeof(*z1[n]) * N);
    }

    /*  Set x and y to the starting value and grab the current time.          */
    x = start;
    y = start;
    t1 = clock();

    /*  Perform the calculation with the C99 library function copysign.       */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z0[m][n] = copysign(x, y);
            y += dy;
        }
        /*  End of y for-loop.                                                */
        x += dx;
    }
    /*  End of x for-loop.                                                    */

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to second, use the macro CLOCKS_PER_SEC    *
     *  provided in time.h.                                                   */
    printf("C99:          %f\n", (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the rss_ringoccs function.                */
    x = start;
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using rssringoccs_Double_Copysign.            */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z1[m][n] = rssringoccs_Double_Copysign(x, y);
            y += dy;
        }
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
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            /*  We'll use the standard library function to check the error.   */
            temp = fabs(z0[m][n] - z1[m][n]);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for loop computing |y0-y1|.                                    */

    /*  Print out the error to 16 decimals (assumes 64-bit precision).        */
    printf("Max Error: %.16f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */

    for (n=0; n<N; ++n)
    {
        free(z0[n]);
        free(z1[n]);
    }

    free(z0);
    free(z1);
    return 0;
}
/*  End of main.                                                              */
