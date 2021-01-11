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

/*  This program plots the real and imaginary parts of the comple exponential *
 *  function using a color gradient to represent the values.                  */

#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

int main(void)
{
    /*  The number of pixels in the x and y axes.                             */
    const unsigned int size = 1024;

    /* Values for the min and max of the x and y axes.                        */
    const double x_min = -4.0;
    const double x_max =  4.0;
    const double y_min = -4.0;
    const double y_max =  4.0;

    /*  Use rssringoccs_Easy_Complex_Plots to produce the plots.              */
    rssringoccs_Easy_Complex_Plots("complex_exp", rssringoccs_Complex_Exp,
                                   size, size, x_min, x_max, y_min, y_max);
    return 0;
}

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_exp_plots.c -o test -lrssringoccs                         *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_src.sh), then change       *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_exp_example.c -o test -lrssringoccs                   *
 *                                                                            *
 *  Running ./test will generate two figure:                                  *
 *      complex_exp_real_part.ppm                                             *
 *      complex_exp_imag_part.ppm                                             *
 *  Which you can view using your favorite image tool. MacOS preview works    *
 *  via open filename.ppm, and gio works on Linux distributions via           *
 *  gio open filename.ppm.                                                    *
 *                                                                            *
 *  This has been tested with several flags to check for strict compliance to *
 *  the C89 standard. The following compiles on gcc and clang:                *
 *      gcc -pedantic -Wall -Wextra -pedantic-errors -std=c89                 *
 *          -ansi -O3 complex_exp_plots.c -o test -lrssringoccs               *
 ******************************************************************************/
