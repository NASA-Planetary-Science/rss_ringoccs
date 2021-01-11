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
 ******************************************************************************/

/*  We'll need stdio for printing the roots of the Bessel J_0 function.       */
#include <stdio.h>

/*  The Newton-Raphson routine is found here.                                 */
#include <rss_ringoccs/include/rss_ringoccs_numerical.h>

/*  The bessel function J_0 is defined here.                                  */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Create an alias for the Bessel function as "f".                           */
double (*f)(double) = &rssringoccs_Bessel_J0_Double;

/*  Numerically compute the derivative.                                       */
static double f_prime(double x)
{
    double h = 1.0e-8;
    return (f(x+h)-f(x-h))/(2*h);
}

int main(void)
{
    /*  Choose a value close to where we think the root is. By looking at the *
     *  graph of J_0(x) we can guess that a root lies between 2 and 3.        */
    double x0 = 2.0;
    double root;
    unsigned int max_iters = 20;

    /*  Use Newton-Raphson to compute the root.                               */
    root = rssringoccs_Newton_Raphson_Double(x0, f, f_prime, max_iters);

    /*  Print the result.                                                     */
    printf("First positive root of J_0: %f\n", root);
    return 0;
}

/******************************************************************************
 *  We can compile this with:                                                 *
 *      gcc newton_rapshon_bessel_function.c -o test -lrssringoccs            *
 *  You may need to add -I/usr/local/include/ or -L/usr/local/lib/ if you do  *
 *  not have these in your path. They should be, by default.                  *
 *  Running the output executable with ./test we get:                         *
 *      First positive root of J_0: 2.404826                                  *
 *  Which is agreement with j_0,1.                                            *
 ******************************************************************************/
