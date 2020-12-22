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
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <complex.h>

/*  Routine for comparing two real valued functions at single precision.      */
void
rssringoccs_Compare_Float_Funcs(const char *f0_name, float (*f0)(float),
                                const char *f1_name, float (*f1)(float),
                                float start, float end, unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    float x, dx;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    float max_err = 0.0F;
    float temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    float *y0, *y1;

    /*  Declare a dummy variable for indexing.                                */
    unsigned long n;

    /*  Declare variables for computing computation time.                     */
    clock_t t1, t2;

    /*  Error check to make sure the user provided valid inputs.              */
    if (start >= end)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_Float_Funcs\n\n"
             "start is greater than or equal to end.\n"
             "Abort computation.\n");
        exit(0);
    }
    else if (N == 0)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_Float_Funcs\n\n"
             "Input sample size is zero. Aborting computation.\n");
        exit(0);
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
}
/*  End of rssringoccs_Compare_Float_Funcs.                                   */

/*  Routine for comparing two real valued functions at double precision.      */
void
rssringoccs_Compare_Double_Funcs(const char *f0_name, double (*f0)(double),
                                 const char *f1_name, double (*f1)(double),
                                 double start, double end, unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    double x, dx;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    double max_err = 0.0;
    double temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    double *y0, *y1;

    /*  Declare a dummy variable for indexing.                                */
    unsigned long n;

    /*  Declare variables for computing computation time.                     */
    clock_t t1, t2;

    /*  Error check to make sure the user provided valid inputs.              */
    if (start >= end)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_Double_Funcs\n\n"
             "start is greater than or equal to end.\n"
             "Abort computation.\n");
        exit(0);
    }
    else if (N == 0)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_Double_Funcs\n\n"
             "Input sample size is zero. Aborting computation.\n");
        exit(0);
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
        temp = fabs(y0[n] - y1[n]);

        /*  Check if the error got larger and set max_err accordingly.        */
        if (max_err < temp)
            max_err = temp;
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 16 decimals (assumes 64-bit precision).        */
    printf("Max Error: %.16f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    free(y0);
    free(y1);
}
/*  End of rssringoccs_Compare_Double_Funcs.                                  */

/*  Routine for comparing two real valued functions at long double precision. */
void
rssringoccs_Compare_LDouble_Funcs(const char *f0_name,
                                  long double (*f0)(long double),
                                  const char *f1_name,
                                  long double (*f1)(long double),
                                  long double start,
                                  long double end,
                                  unsigned long N)
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
        exit(0);
    }
    else if (N == 0)
    {
        puts("\nError Encountered: rss_ringoccs\n"
             "\r\trssringoccs_Compare_LDouble_Funcs\n\n"
             "Input sample size is zero. Aborting computation.\n");
        exit(0);
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

/*  Routine for comparing two complex valued functions at single precision.   */
void
rssringoccs_Compare_CFloat_Funcs(const char *f0_name,
                                 rssringoccs_ComplexFloat
                                   (*f0)(rssringoccs_ComplexFloat),
                                 const char *f1_name,
                                 float _Complex
                                   (*f1)(float _Complex),
                                 float start,
                                 float end,
                                 unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    float x, y, ds;
    float x_c, y_c, x_s, y_s, dx, dy;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    float max_err = 0.0F;
    float temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
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

    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z0[m][n] = f0(rssringoccs_CFloat_Rect(x, y));
            x += ds;
        }
        y += ds;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    x = start;
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z1[m][n] = f1(x + _Complex_I*y);
            x += ds;
        }
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
            x_s = rssringoccs_CFloat_Real_Part(z0[m][n]);
            y_s = rssringoccs_CFloat_Imag_Part(z0[m][n]);

            x_c = crealf(z1[m][n]);
            y_c = cimagf(z1[m][n]);

            dx = x_s - x_c;
            dy = y_s - y_c;

            /*  We'll use the standard library function to check the error.   */
            temp = sqrtf(dx*dx + dy*dy);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.8f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_CFloat_Funcs.                                  */

void
rssringoccs_Compare_CDouble_Funcs(const char *f0_name,
                                  rssringoccs_ComplexDouble
                                    (*f0)(rssringoccs_ComplexDouble),
                                  const char *f1_name,
                                  double _Complex
                                    (*f1)(double _Complex),
                                  double start, double end,
                                  unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    double x, y, ds;
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

    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z0[m][n] = f0(rssringoccs_CDouble_Rect(x, y));
            x += ds;
        }
        y += ds;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    x = start;
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z1[m][n] = f1(x + (double _Complex)_Complex_I * y);
            x += ds;
        }
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
            x_s = rssringoccs_CDouble_Real_Part(z0[m][n]);
            y_s = rssringoccs_CDouble_Imag_Part(z0[m][n]);

            x_c = creal(z1[m][n]);
            y_c = cimag(z1[m][n]);

            dx = x_s - x_c;
            dy = y_s - y_c;

            /*  We'll use the standard library function to check the error.   */
            temp = sqrt(dx*dx + dy*dy);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.16f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_CDouble_Funcs.                                 */

extern void
rssringoccs_Compare_CLDouble_Funcs(const char *f0_name,
                                   rssringoccs_ComplexLongDouble
                                     (*f0)(rssringoccs_ComplexLongDouble),
                                   const char *f1_name,
                                   long double _Complex
                                     (*f1)(long double _Complex),
                                   long double start,
                                   long double end,
                                   unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    long double x, y, ds;
    long double x_c, y_c, x_s, y_s, dx, dy;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    long double max_err = 0.0L;
    long double temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    rssringoccs_ComplexLongDouble **z0;
    long double _Complex **z1;

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

    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z0[m][n] = f0(rssringoccs_CLDouble_Rect(x, y));
            x += ds;
        }
        y += ds;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    x = start;
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z1[m][n] = f1(x + (long double _Complex)_Complex_I*y);
            x += ds;
        }
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
            x_s = rssringoccs_CLDouble_Real_Part(z0[m][n]);
            y_s = rssringoccs_CLDouble_Imag_Part(z0[m][n]);

            x_c = creall(z1[m][n]);
            y_c = cimagl(z1[m][n]);

            dx = x_s - x_c;
            dy = y_s - y_c;

            /*  We'll use the standard library function to check the error.   */
            temp = sqrtl(dx*dx + dy*dy);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.24Lf\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_CLDouble_Funcs.                                */

void
rssringoccs_Compare_Real_CFloat_Funcs(const char *f0_name,
                                      float (*f0)(rssringoccs_ComplexFloat),
                                      const char *f1_name,
                                      float (*f1)(float _Complex),
                                      float start, float end, unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    float x, y, ds;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    float max_err = 0.0F;
    float temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    float **z0;
    float **z1;

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

    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z0[m][n] = f0(rssringoccs_CFloat_Rect(x, y));
            x += ds;
        }
        y += ds;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    x = start;
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z1[m][n] = f1(x + _Complex_I*y);
            x += ds;
        }
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
            temp = fabsf(z0[m][n] - z1[m][n]);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.8f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_Real_CFloat_Funcs.                             */

void
rssringoccs_Compare_Real_CDouble_Funcs(const char *f0_name,
                                       double (*f0)(rssringoccs_ComplexDouble),
                                       const char *f1_name,
                                       double (*f1)(double _Complex),
                                       double start, double end,
                                       unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    double x, y, ds;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    double max_err = 0.0;
    double temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    double **z0;
    double **z1;

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

    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z0[m][n] = f0(rssringoccs_CDouble_Rect(x, y));
            x += ds;
        }
        y += ds;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    x = start;
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z1[m][n] = f1(x + (double)_Complex_I*y);
            x += ds;
        }
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
            temp = fabs(z0[m][n] - z1[m][n]);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.16f\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_Real_CDouble_Funcs.                            */

void
rssringoccs_Compare_Real_CLDouble_Funcs(const char *f0_name,
                                        long double
                                          (*f0)(rssringoccs_ComplexLongDouble),
                                        const char *f1_name,
                                        long double (*f1)(long double _Complex),
                                        long double start, long double end,
                                        unsigned long N)
{
    /*  Declare variables for sampling the region [start, end].               */
    long double x, y, ds;

    /*  Declare variables for computing the maximum difference between the    *
     *  two provided functions f0 and f1.                                     */
    long double max_err = 0.0L;
    long double temp;

    /*  Declare two pointers to represent arrays for f0(x) and f1(x).         */
    long double **z0;
    long double **z1;

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

    for (m=0; m<N; ++m)
    {
        z0[m] = malloc(sizeof(*z0[m]) * N);
        z1[m] = malloc(sizeof(*z1[m]) * N);
    }

    /*  Set x to the starting value and grab the current time.                */
    x = start;
    y = start;
    t1 = clock();

    /*  Perform the calculation for the f0 function.                          */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z0[m][n] = f0(rssringoccs_CLDouble_Rect(x, y));
            x += ds;
        }
        y += ds;
    }

    /*  Grab the current clock time again.                                    */
    t2 = clock();

    /*  t2-t1 is the number of clock cycles that have passed between grabbing *
     *  t1 and t2. To convert this to seconds, use the macro CLOCKS_PER_SEC   *
     *  provided in time.h.                                                   */
    printf("%s: %f\n", f0_name, (double)(t2-t1)/CLOCKS_PER_SEC);

    /*  Restart the computation for the f1 function.                          */
    x = start;
    y = start;

    /*  Reset the clock.                                                      */
    t1 = clock();

    /*  Perform the computation using f1 instead of f0.                       */
    for (m=0; m<N; ++m)
    {
        for (n=0; n<N; ++n)
        {
            z1[m][n] = f1(x + (long double)_Complex_I*y);
            x += ds;
        }
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
            temp = fabsl(z0[m][n] - z1[m][n]);

            /*  Check if the error got larger and set max_err accordingly.    */
            if (max_err < temp)
                max_err = temp;
        }
    }
    /*  End of for-loop computing |y0-y1|.                                    */

    /*  Print out the error to 8 decimals (assumes 32-bit precision).         */
    printf("Max Error: %.24Lf\n", max_err);

    /*  Free the pointers we've malloc'd.                                     */
    for (m=0; m<N; ++m)
    {
        free(z0[m]);
        free(z1[m]);
    }

    free(z0);
    free(z1);
}
/*  End of rssringoccs_Compare_Real_CLDouble_Funcs.                           */
