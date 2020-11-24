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

/*  We'll need the cosh function.                                             */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  If rss_ringoccs built correctly, rss_ringoccs_complex.h is located in     *
 *  /usr/local/include/rss_ringoccs/include. We can include this as follows:  */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  The number of pixels in the x and y axes.                                 */
const unsigned int size = 1024;

/*  This function is used to set the current pixel of the output ppm to the   *
 *  desired color. fp is the filename we'll be using.                         */
static void color(double mag, FILE *fp)
{
    unsigned char red, green, blue;

    /*  Use an RGB rainbow gradient to color the current pixel.               */
    if (mag < 64)
    {
        red   = (unsigned char)0;
        green = (unsigned char)4*mag;
        blue  = (unsigned char)255;
    }
    else if (mag < 128)
    {
        red   = (unsigned char)0;
        green = (unsigned char)255;
        blue  = (unsigned char)(255 - 4*(mag - 64));
    }
    else if (mag < 192)
    {
        red   = (unsigned char)4*(mag-128);
        green = (unsigned char)255;
        blue  = (unsigned char)0;
    }
    else
    {
        red   = (unsigned char)255;
        green = (unsigned char)(255 - 4*(mag-192));
        blue  = (unsigned char)0;
    }

    /*  Color the current pixel.                                              */
    fputc(red,   fp);
    fputc(green, fp);
    fputc(blue,  fp);
}

int main(void)
{
    /* Values for the min and max of the x and y axes.                        */
    double x_min = -4.0;
    double x_max =  4.0;
    double y_min = -4.0;
    double y_max =  4.0;

    /*  Variables for the min and max of the real and imaginary parts of      *
     *  cos(z) in the region of interest.                                     */
    double min_re, max_re, min_im, max_im;

    /*  Declare variables needed for the computation.                         */
    unsigned int x, y;
    double z_x, z_y, w_x, w_y;
    rssringoccs_ComplexDouble z, w;

    /*  Declare a variable for the output files.                              */
    FILE **fp;

    /*  There are two files, so malloc memory for two.                        */
    fp = malloc(sizeof(*fp) * 2);

    /*  Create the files and give them write permissions.                     */
    fp[0] = fopen("complex_cos_real_part.ppm", "w");
    fp[1] = fopen("complex_cos_imag_part.ppm", "w");

    /*  Find the min and max of the function so we can set the color scale.   */
    if (fabs(y_min) < fabs(y_max))
    {
        max_re = rssringoccs_Cosh_Double(y_max);
        min_re = -max_re;

        max_im = rssringoccs_Cosh_Double(fabs(y_max));
        min_im = -max_re;
    }
    else
    {
        max_re = rssringoccs_Cosh_Double(y_min);
        min_re = -max_re;
        max_im = rssringoccs_Cosh_Double(fabs(y_min));
        min_im = -max_re;
    }

    /*  Needed to create the output ppm file. This is the preamble.           */
    fprintf(fp[0], "P6\n%d %d\n255\n", size, size);
    fprintf(fp[1], "P6\n%d %d\n255\n", size, size);

    /*  Compress the values to a smaller range to bring out some dynamics in  *
     *  the image.                                                            */
    min_re = atan(min_re);
    max_re = atan(max_re);
    min_im = atan(min_im);
    max_im = atan(max_im);

    /*  Translate max_im and max_re to be positive.                           */
    max_re = max_re - min_re;
    max_im = max_im - min_im;

    /*  Loop over each pixel and color it based on the value of cos(x+iy).    */
    for (y=0; y<size; ++y)
    {
        /*  We want to center z_y so scale and shift. This makes the output   *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        z_y = y * (y_max - y_min)/(size - 1) + y_min;

        for (x=0; x<size; ++x)
        {
            /*  Similarly, center z_x.                                        */
            z_x = x * (x_max - x_min)/(size - 1) + x_min;

            /*  Set z to x+iy.                                                */
            z = rssringoccs_Complex_Rect(z_x, z_y);

            /*  Compute the complex cosine of z.                              */
            w = rssringoccs_Complex_Cos(z);

            /*  Extract the real and imaginary parts of w.                    */
            w_x = rssringoccs_Complex_Real_Part(w);
            w_y = rssringoccs_Complex_Imag_Part(w);


            /*  Scale w_x and w_y to lie in the range [0, 255].               */
            w_x = 255*(atan(w_x) - min_re)/max_re;
            w_y = 255*(atan(w_y) - min_im)/max_im;

            /*  Color the current pixel.                                      */
            color(w_x, fp[0]);
            color(w_y, fp[1]);
        }
    }

    /*  Free the memory allocated to fp.                                      */
    free(fp);
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
