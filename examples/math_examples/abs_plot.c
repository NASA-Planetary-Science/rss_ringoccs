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
 *      Plot the absolute value function y = |x|.                             *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 11, 2020                                             *
 ******************************************************************************/

/*  We'll use these macros for drawing the figure in a pbm file later.        */
#define BLACK (unsigned char)0
#define WHITE (unsigned char)255

/*  Needed for creating the output file.                                      */
#include <stdio.h>

/*  Contains malloc and free.                                                 */
#include <stdlib.h>

/*  The absolute value function is found here.                                */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  The number of pixels in the x and y axes.                                 */
const unsigned int size = 1024;

int main(void)
{
    /* Values for the min and max of the x and y axes.                        */
    double x_min =  -1.2;
    double x_max =  1.2;
    double y_min = -0.2;
    double y_max =  1.2;

    /*  Set a parameter for the thickness of the curve and the axes.          */
    double pixel_width = 0.002;

    /*  Declare variables needed for the computation.                         */
    unsigned int x, y;
    double Px, Py, *abs_x, rcp_factor, diff, abs_Px, abs_Py;

    /*  Declare a variable for the output file.                               */
    FILE *fp;

    /*  Silly check to make sure the user provided a valid range for x and y. */
    if (x_max <= x_min)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tabs_plot.c\n\n"
             "x_min is greater than or equal to x_max.\n");
        exit(0);
    }
    else if (y_max <= y_min)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tabs_plot.c\n\n"
             "y_min is greater than or equal to y_max.\n");
        exit(0);
    }

    /*  Another silly error check to make sure size is greater than 1.        */
    if (size == 0)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tabs_plot.c\n\n"
             "Input size is zero. Aborting computation.\n");
        exit(0);
    }
    else if (size == 1)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tabs_plot.c\n\n"
             "Input size is one. This will cause divide-by-zero.\n"
             "Aborting computation.\n");
        exit(0);
    }

    /*  Otherwise, set rcp_factor to 1/(size-1) so that we can scale between  *
     *  the pixel grid [0, size] x [0, size] and the Cartesian coordinates    *
     *  [x_min, x_max] x [y_min, y_max].                                      */
    else
        rcp_factor = 1.0 / (size - 1.0);

    /*  Create the files and give them write permissions.                     */
    fp = fopen("absolute_value_plot.pgm", "w");

    /*  Needed to create the output pbm file. This is the preamble.           */
    fprintf(fp, "P5\n%d %d\n255\n", size, size);

    /*  Allocate memory for the two variables.                                */
    abs_x = malloc(sizeof(*abs_x) * size);

    /*  Loop through and compute abs_arg = |arg|.                             */
    for (x=0; x<size; ++x)
    {
        Px = x * (x_max - x_min) * rcp_factor + x_min;
        abs_x[x] = rssringoccs_Double_Abs(Px);
    }

    /*  Loop over each pixel and color it based on the value of cos(x+iy).    */
    for (y=0; y<size; ++y)
    {
        /*  We want to center Py so scale and shift. This makes the output    *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        Py = (size-y) * (y_max - y_min) * rcp_factor + y_min;
        abs_Py = rssringoccs_Double_Abs(Py);

        for (x=0; x<size; ++x)
        {
            /*  Similarly, center Px.                                         */
            Px = x * (x_max - x_min) * rcp_factor + x_min;
            diff = rssringoccs_Double_Abs(Py - abs_x[x]);
            abs_Px = rssringoccs_Double_Abs(Px);

            /*  Color in the x-axis.                                          */
            if (abs_Px < pixel_width)
                fputc(WHITE, fp);

            /*  Color in the y-axis.                                          */
            else if (abs_Py < pixel_width)
                fputc(WHITE, fp);

            /*  Color in the function y = |x|.                                */
            else if (diff < pixel_width)
                fputc(WHITE, fp);

            /*  Otherwise, color the background black.                        */
            else
                fputc(BLACK, fp);
        }
    }
    return 0;
}

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc abs_plot.c -o test -lrssringoccs                                  *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_librssringoccs.sh), then   *
 *  change the -L option to the correct location. If /usr/local/include/ is   *
 *  not in your path, add the -I option as follows:                           *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              abs_plot.c -o test -lrssringoccs                              *
 *                                                                            *
 *  Running ./test will generate the figure absolute_value_plot.pgm which you *
 *  can view using your favorite image tool. MacOS preview works via open     *
 *  filename.ppm, and gio works on Linux distributions via gio open           *
 *  filename.ppm.                                                             *
 *                                                                            *
 *  This has been tested with several flags to check for strict compliance to *
 *  the C89 standard. The following compiles on gcc and clang:                *
 *      gcc -pedantic -Wall -Wextra -pedantic-errors -std=c89 -Wconversion    *
 *          -ansi -O3 abs_plot.c -o test -lrssringoccs                        *
 ******************************************************************************/
