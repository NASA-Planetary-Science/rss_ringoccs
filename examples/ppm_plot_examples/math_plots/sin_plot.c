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

/*  Needed for creating the output file.                                      */
#include <stdio.h>

/*  Contains malloc and free.                                                 */
#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  The number of pixels in the x and y axes.                                 */
const unsigned int size = 1024;

int main(void)
{
    /* Values for the min and max of the x and y axes.                        */
    double x_min =  -0.2;
    double x_max =  rssringoccs_Two_Pi;
    double y_min = -1.2;
    double y_max =  1.2;

    double y_to_x_scale = (x_max - x_min)/(y_max - y_min);
    double pixel_width = 0.002;

    /*  Declare variables needed for the computation.                         */
    unsigned int x, y;
    double Px, Py, erfcx_Px, rcp_factor;

    /*  Declare a variable for the output files.                              */
    FILE *fp;

    /*  Silly check to make sure the user provided a valid range for x and y. */
    if (x_max <= x_min)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tcomplex_cos_plots.c\n\n"
             "x_min is greater than or equal to x_max.\n");
        exit(0);
    }

    /*  Another silly error check to make sure size is greater than 1.        */
    if (size == 0)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tcomplex_cos_plots.c\n\n"
             "Input size is zero. Aborting computation.\n");
        exit(0);
    }
    else if (size == 1)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\tcomplex_cos_plots.c\n\n"
             "Input size is one. This will cause divide-by-zero.\n"
             "Aborting computation.\n");
        exit(0);
    }

    /*  Create the files and give them write permissions.                     */
    fp = fopen("erfcx_plot.pgm", "w");

    /*  Needed to create the output ppm file. This is the preamble.           */
    fprintf(fp, "P5\n%d %d\n255\n", size, size);

    /*  To translate from the pixel (x, y) to the point (z_x, z_y) lying in   *
     *  the rectangle [x_min, x_max] x [y_min, y_max] we'll need this term.   */
    rcp_factor = 1.0 / (size - 1.0);

    /*  Loop over each pixel and color it based on the value of cos(x+iy).    */
    for (y=0; y<size; ++y)
    {
        /*  We want to center Px so scale and shift. This makes the output    *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        Py = (size-y) * (y_max - y_min) * rcp_factor + y_min;

        for (x=0; x<size; ++x)
        {
            /*  Similarly, center Py.                                         */
            Px = x * (x_max - x_min) * rcp_factor + x_min;
            erfcx_Px = rssringoccs_Double_Sin(Px);

            if ((rssringoccs_Double_Abs(Px) < pixel_width*y_to_x_scale) ||
                (rssringoccs_Double_Abs(Py - erfcx_Px) < pixel_width) ||
                (rssringoccs_Double_Abs(Py) < pixel_width))
                fputc((unsigned char)255, fp);
            else
                fputc((unsigned char)0, fp);
        }
    }
    return 0;
}

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc complex_cos_plots.c -o test -lrssringoccs                         *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_src.sh), then change       *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_cos_example.c -o test -lrssringoccs                   *
 *                                                                            *
 *  Running ./test will generate two figure:                                  *
 *      complex_cos_real_part.ppm                                             *
 *      complex_cos_imag_part.ppm                                             *
 *  Which you can view using your favorite image tool. MacOS preview works    *
 *  via open filename.ppm, and gio works on Linux distributions via           *
 *  gio open filename.ppm.                                                    *
 *                                                                            *
 *  This has been tested with several flags to check for strict compliance to *
 *  the C89 standard. The following compiles on gcc and clang:                *
 *      gcc -pedantic -Wall -Wextra -pedantic-errors -std=c89                 *
 *          -ansi -O3 complex_cos_plots.c -o test -lrssringoccs               *
 ******************************************************************************/
