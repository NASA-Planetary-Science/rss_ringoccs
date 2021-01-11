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

/*  Routine for comparing two complex valued functions at single precision.   */
void
rssringoccs_Compare_CFloat_Funcs(
    const char *f0_name,
    rssringoccs_ComplexFloat (*f0)(rssringoccs_ComplexFloat),
    const char *f1_name,
    float _Complex (*f1)(float _Complex),
    const float start,
    const float end,
    const unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    float x, y, ds;

    /*  And declare variables for computing the |z0 - z1|.                    */
    float x_c, y_c, x_s, y_s, dx, dy;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    float max_err = 0.0F;
    float temp;

    /*  Declare two pointers to represent arrays for f0(x+iy) and f1(x+iy).   */
    rssringoccs_ComplexFloat **z0;
    float _Complex **z1;

    /*  Declare a dummy variable for indexing.                                */
    unsigned long m, n;

    /*  Declare variables for computing computation time.                     */
    clock_t t1, t2;

    /*  Error check to make sure the user provided valid inputs.              */
    if (start >= end)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_CFloat_Funcs\n\n"
             "start is greater than or equal to end.\n"
             "Abort computation.\n");
        exit(0);
    }
    else if (N == 0)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_CFloat_Funcs\n\n"
             "Input sample size is zero. Aborting computation.\n");
        exit(0);
    }

    /*  We'll increment evenly throughout the region.                         */
    ds = (end - start) / N;

    /*  Allocate memory for the two pointers we've declared.                  */
    z0 = malloc(sizeof(*z0) * N);
    z1 = malloc(sizeof(*z1) * N);

    /*  We've allocated memory for a pointer to N pointers for z0 and z1. Set *
     *  each element to a block of memory of size N. This will mimic a matrix *
     *  or array of NxN complex numbers. We can access the (m,n) entry via    *
     *  z0[m][n] and z1[m][n].                                                */
    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }
    /*  Be sure to free all of this when we're done!                          */

    /*  Set y to the starting value and grab the current time.                */
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        /*  Set x to the start value and loop over the next roow.             */
        x = start;
        for (n=0; n<N; ++n)
        {
            z0[m][n] = f0(rssringoccs_CFloat_Rect(x, y));

            /*  Move to the next row entry in the same column.                */
            x += ds;
        }

        /*  Increment y to the next column.                                   */
        y += ds;
    }
    /*  End of for loop computing z0 = f0(x+iy).                              */

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
        /*  Reset x to the start of the row.                                  */
        x = start;
        for (n=0; n<N; ++n)
        {
            z1[m][n] = f1(x + _Complex_I*y);
            x += ds;
        }

        /*  Increment y to the next column.                                   */
        y += ds;
    }
    /*  End of for loop computing z1 = f1(x+iy).                              */

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
            /*  Extract the real and imaginary parts from z0.                 */
            x_s = rssringoccs_CFloat_Real_Part(z0[m][n]);
            y_s = rssringoccs_CFloat_Imag_Part(z0[m][n]);

            /*  Extract the real and imaginary parts from z1.                 */
            x_c = crealf(z1[m][n]);
            y_c = cimagf(z1[m][n]);

            /*  We'll compute |z0-z1| via the Pythagorean formula. Compute    *
             *  the distance between the z0 and z1 x components and y         *
             *  components, respectively.                                     */
            dx = x_s - x_c;
            dy = y_s - y_c;

            /*  We'll use the standard library function to check the error.   *
             *  Note, this is the absolute error, not relative.               */
            temp = sqrtf(dx*dx + dy*dy);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |z0-z1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.8f\n", (double)max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    /*  And free the pointers to the pointers (z0 and z1).                    */
    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_CFloat_Funcs.                                  */
