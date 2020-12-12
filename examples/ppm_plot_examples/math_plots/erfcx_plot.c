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
 *      Plot the scaled complementary error function erfcx(x)                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  The absolute value function is found here.                                */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Plotting routines defined here.                                           */
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

/*  Routine for plotting the scaled complementary error function.             */
int main(void)
{
    /* Values for the min and max of the x and y axes.                        */
    double x_min = -0.6;
    double x_max =  4.2;
    double y_min = -0.6;
    double y_max =  1.8;

    /*  The number of pixels in the x axes.                                   */
    const unsigned int x_size = 4*1024;

    /*  The number of pixels in the y axes.                                   */
    const unsigned int y_size = 2*1024;

    /*  Plot the figure using the rss_ringoccs ppm_plot routines.             */
    rssringoccs_Easy_Real_Plots("erfcx", rssringoccs_Double_Erfcx,
                                 x_size, y_size, x_min, x_max, y_min, y_max);

    return 0;
}
/*  End of main.                                                              */

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc erfcx_plot.c -o test -lrssringoccs                                *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_src.sh), then change       *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              erfcx_plot.c -o test -lrssringoccs                            *
 *                                                                            *
 *  Running ./test will generate the figure erfcx_plot.pgm which you can view *
 *  using your favorite image tool. MacOS preview works via open filename.pgm *
 *  and gio works on Linux distributions via gio open filename.pgm.           *
 *                                                                            *
 *  This has been tested with several flags to check for strict compliance to *
 *  the C89 standard. The following compiles on gcc and clang:                *
 *      gcc -pedantic -Wall -Wextra -pedantic-errors -std=c89                 *
 *          -ansi -O3 erfcx_plot.c -o test -lrssringoccs                      *
 ******************************************************************************/
