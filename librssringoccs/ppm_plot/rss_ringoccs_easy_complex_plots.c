#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

void rssringoccs_Easy_Complex_Plots(
    const char *func_name,
    rssringoccs_ComplexDouble(*func)(rssringoccs_ComplexDouble),
    unsigned int x_size, unsigned int y_size,
    const double x_min, const double x_max,
    const double y_min, const double y_max
)
{
    unsigned int x, y;
    size_t string_length;
    FILE **fp;
    char *re_file, *im_file, *arg_file, *mag_file;
    double x_rcp_factor, y_rcp_factor;
    rssringoccs_ComplexDouble z, w;
    double z_x, z_y, w_x, w_y;
    double arg, mag;

    if (func_name == NULL)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Easy_Complex_Plots\n\n"
             "Input func_name string is NULL. Aborting.\n");
        exit(0);
    }

    /*  Silly check to make sure the user provided a valid range for x and y. */
    if (x_max <= x_min)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Easy_Complex_Plots.c\n\n"
             "x_min is greater than or equal to x_max.\n");
        exit(0);
    }
    else if (y_max <= y_min)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Easy_Complex_Plots.c\n\n"
             "y_min is greater than or equal to y_max.\n");
        exit(0);
    }

    /*  Another silly error check to make sure size is greater than 1.        */
    if ((x_size == 0) || (y_size == 0))
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Easy_Complex_Plots.c\n\n"
             "Input size is zero. Aborting computation.\n");
        exit(0);
    }
    else if ((x_size == 1) || (y_size == 1))
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\trssringoccs_Easy_Complex_Plots.c\n\n"
             "Input size is one. This will cause divide-by-zero.\n"
             "Aborting computation.\n");
        exit(0);
    }
    string_length = strlen(func_name);
    re_file = malloc(sizeof(*re_file) *(string_length + 10));
    im_file = malloc(sizeof(*im_file) *(string_length + 10));
    mag_file = malloc(sizeof(*re_file) *(string_length + 9));
    arg_file = malloc(sizeof(*im_file) *(string_length + 9));

    x_rcp_factor = 1.0/(x_size - 1.0);
    y_rcp_factor = 1.0/(y_size - 1.0);

    strcpy(re_file, func_name);
    strcpy(im_file, func_name);
    strcpy(mag_file, func_name);
    strcpy(arg_file, func_name);
    strcat(re_file, "_real.ppm");
    strcat(im_file, "_imag.ppm");
    strcat(mag_file, "_mag.ppm");
    strcat(arg_file, "_arg.ppm");

    fp = malloc(sizeof(*fp) * 2);

    fp[0] = fopen(re_file, "w");
    fp[1] = fopen(im_file, "w");
    fp[3] = fopen(mag_file, "w");
    fp[4] = fopen(arg_file, "w");

    /*  Needed to create the output ppm file. This is the preamble.           */
    fprintf(fp[0], "P6\n%d %d\n255\n", x_size, y_size);
    fprintf(fp[1], "P6\n%d %d\n255\n", x_size, y_size);
    fprintf(fp[3], "P6\n%d %d\n255\n", x_size, y_size);
    fprintf(fp[4], "P6\n%d %d\n255\n", x_size, y_size);

    for (y=0; y<y_size; ++y)
    {
        /*  We want to center z_y so scale and shift. This makes the output   *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        z_y = (y_size - y - 1.0) * (y_max - y_min) * y_rcp_factor + y_min;

        for (x=0; x<x_size; ++x)
        {
            /*  Similarly, center z_x.                                        */
            z_x = x * (x_max - x_min) * x_rcp_factor + x_min;

            /*  Set z to x+iy.                                                */
            z = rssringoccs_Complex_Rect(z_x, z_y);

            /*  Compute the complex cosine of z.                              */
            w = func(z);

            /*  Extract the real and imaginary parts of w.                    */
            w_x = rssringoccs_Complex_Real_Part(w);
            w_y = rssringoccs_Complex_Imag_Part(w);
            arg = rssringoccs_Complex_Argument(w);
            mag = rssringoccs_Complex_Abs(w);

            /*  Color the current pixel.                                      */
            rssringoccs_RGB_Scaled_Gradient(w_x, fp[0]);
            rssringoccs_RGB_Scaled_Gradient(w_y, fp[1]);
            rssringoccs_RGB_Scaled_Gradient(mag, fp[3]);
            rssringoccs_RGB_Linear_Gradient(arg, -ONE_PI, ONE_PI, fp[4]);
        }
    }

    free(fp);
    free(re_file);
    free(im_file);
    free(mag_file);
    free(arg_file);
}
