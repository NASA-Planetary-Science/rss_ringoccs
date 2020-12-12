#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

/*  We'll use these macros for drawing the figure in a pgm file later.        */
#define BLACK (unsigned char)0
#define WHITE (unsigned char)255

/*  Routine for plotting the absolute value function.                         */
void rssringoccs_Easy_Real_Plots(const char *func_name, double (*f)(double),
                                 unsigned int x_size, unsigned int y_size,
                                 const double x_min, const double x_max,
                                 const double y_min, const double y_max)
{

    /*  Set a parameter for the thickness of the curve and the axes.          */
    double pixel_width = 0.002;

    /*  Declare variables for the pixel (x, y).                               */
    unsigned int x, y;

    /*  Declare variables for parsing the input function name and creating    *
     *  the output file name.                                                 */
    size_t string_length;
    char *filename;

    /*  Declare variables for the Cartesian coordinates corresponding to the  *
     *  pixel (x, y). These will represent arrays [x_min, x_max] and          *
     *  [y_min, y_max], respectively.                                         */
    double *Px, *Py;

    /*  Declare variables for |Px| and |Py|, respectively.                    */
    double *abs_Px, *abs_Py;

    /*  Declare a pointer for f(x) where x varies over [x_min, x_max].        */
    double *f_of_x;

    /*  Lastly, declare two variables used later in the computation.          */
    double rcp_factor_x, rcp_factor_y, diff;

    /*  Declare a variable for the output file.                               */
    FILE *fp;

    /*  Silly check to make sure the user provided a valid range for x and y. */
    if (x_max <= x_min)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\r\tssringoccs_Easy_Real_Plots\n\n"
             "x_min is greater than or equal to x_max.\n");
        exit(0);
    }
    else if (y_max <= y_min)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\r\tssringoccs_Easy_Real_Plots\n\n"
             "y_min is greater than or equal to y_max.\n");
        exit(0);
    }

    /*  Another silly error check to make sure size is greater than 1.        */
    else if (x_size == 0)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\r\tssringoccs_Easy_Real_Plots\n\n"
             "Input x size is zero. Aborting computation.\n");
        exit(0);
    }
    else if (y_size == 0)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\r\tssringoccs_Easy_Real_Plots\n\n"
             "Input y size is zero. Aborting computation.\n");
        exit(0);
    }
    else if (x_size == 1)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\r\tssringoccs_Easy_Real_Plots\n\n"
             "Input x size is one. This will cause divide-by-zero.\n"
             "Aborting computation.\n");
        exit(0);
    }
    else if (y_size == 1)
    {
        puts("Error Encountered: rss_ringoccs\n"
             "\r\tssringoccs_Easy_Real_Plots\n\n"
             "Input y size is one. This will cause divide-by-zero.\n"
             "Aborting computation.\n");
        exit(0);
    }

    /*  Otherwise, set rcp_factor to 1/(size-1) so that we can scale between  *
     *  the pixel grid [0, size] x [0, size] and the Cartesian coordinates    *
     *  [x_min, x_max] x [y_min, y_max].                                      */
    rcp_factor_x = 1.0 / (x_size - 1.0);
    rcp_factor_y = 1.0 / (y_size - 1.0);

    /*  Get the length of the input string func_name. strlen is defined in    *
     *  the standard library string.h.                                        */
    string_length = strlen(func_name);

    /*  Allocate enough memory for filename to include func_name plus the     *
     *  suffix "_plot.pgm" plus the NULL terminator at the end. So we need    *
     *  string_length + 9 + 1 in total.                                       */
    filename = malloc(sizeof(*filename) *(string_length + 10));

    /*  Copy func_name to filename and then concatenate "_plot.pgm" to it.    */
    strcpy(filename, func_name);
    strcat(filename, "_plot.pgm");

    /*  Create the file and give it write permissions.                        */
    fp = fopen(filename, "w");

    /*  Needed to create the output pgm file. This is the preamble.           */
    fprintf(fp, "P5\n%d %d\n255\n", x_size, y_size);

    /*  Allocate memory for the five variables.                               */
    Px     = malloc(sizeof(Px) * x_size);
    Py     = malloc(sizeof(Py) * y_size);
    abs_Px = malloc(sizeof(*abs_Px) * x_size);
    abs_Py = malloc(sizeof(*abs_Py) * y_size);
    f_of_x = malloc(sizeof(*f_of_x) * x_size);

    /*  Loop through and compute the Cartesian x coordinate corresponding to  *
     *  a given integer in the range [0, x_size]. Compute the absolute value  *
     *  of this as well.                                                      */
    for (x=0; x<x_size; ++x)
    {
        /*  We want to center Px so scale and shift. This makes the output    *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        Px[x] = x * (x_max - x_min) * rcp_factor_x + x_min;

        /*  Compute the absolute value for plotting the x axis later.         */
        abs_Px[x] = rssringoccs_Double_Abs(Px[x]);
        f_of_x[x] = f(Px[x]);
    }
    /*  End of for loop calculating Px and |Px|.                              */

    /*  Do the same thing for y.                                              */
    for (y=0; y<y_size; ++y)
    {
        /*  Similarly, center Py.                                             */
        Py[y] = (y_size - y - 1.0) * (y_max - y_min) * rcp_factor_y + y_min;

        /*  Compute |Py| to plot the y-axis later.                            */
        abs_Py[y] = rssringoccs_Double_Abs(Py[y]);
    }
    /*  End of for loop calculating Py and |Py|.                              */

    /*  Loop over each pixel and color it based on the value |y - |x||.       */
    for (y=0; y<y_size; ++y)
    {
        /*  Loop over all of the x pixels.                                    */
        for (x=0; x<x_size; ++x)
        {
            /*  Compute the distance between the current y pixel and the      *
             *  pixel corresponding to y = |x|.                               */
            diff = rssringoccs_Double_Abs(Py[y] - f_of_x[x]);

            /*  Color in the x-axis.                                          */
            if (abs_Px[x] < pixel_width)
                fputc(WHITE, fp);

            /*  Color in the y-axis.                                          */
            else if (abs_Py[y] < pixel_width)
                fputc(WHITE, fp);

            /*  Color in the function y = |x|.                                */
            else if (diff < pixel_width)
                fputc(WHITE, fp);

            /*  Otherwise, color the background black.                        */
            else
                fputc(BLACK, fp);
        }
        /*  End of looping over the x variable.                               */
    }
    /*  End of looping over the y variable.                                   */

    /*  Close the file.                                                       */
    fclose(fp);

    /*  Free all of the memory we malloc'd.                                   */
    free(Px);
    free(Py);
    free(abs_Px);
    free(abs_Py);
    free(f_of_x);
    free(filename);
}
/*  End of rssringoccs_Easy_Real_Plots.                                       */
