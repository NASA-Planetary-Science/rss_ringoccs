#include <stdio.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

static rssringoccs_ComplexDouble
square(rssringoccs_ComplexDouble z, rssringoccs_ComplexDouble c)
{
    rssringoccs_ComplexDouble out;

    out = rssringoccs_CDouble_Multiply(z, z);
    out = rssringoccs_CDouble_Add(out, c);
    return out;
}

int main(void)
{
    /*  Declare a variable for the output file and give it write permission.  */
    FILE *fp;
    unsigned int x, y, n;
    unsigned char red, green, blue;
    double z_x, z_y, norm;
    double rcp_factor;
    unsigned int maxIterations = 256;
    rssringoccs_ComplexDouble z, c;

    unsigned int size = 1024;

    const double x_min = -2.0;
    const double x_max =  2.0;
    const double y_min = -2.0;
    const double y_max =  2.0;

    double radius = 4.0;

    fp = fopen("mandelbrot_set.ppm", "w");
    fprintf(fp, "P6\n%d %d\n255\n", size, size);

    rcp_factor = 1.0/(size-1.0);

    /*  Loop through each pixel.                                              */
    for (y=0; y<size; ++y)
    {
        /*  We want to center z_y so scale and shift. This makes the output   *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        z_y = (size - y - 1.0) * (y_max - y_min) * rcp_factor + y_min;

        for (x=0; x<
        size; ++x)
        {
            /*  Similarly, center z_x.                                        */
            z_x = x * (x_max - x_min) * rcp_factor + x_min;

            /*  Compute the complex number z_x + i z_y.                       */
            c = rssringoccs_CDouble_Rect(z_x, z_y);

            /*  Reset starting Real and Imaginary parts to zero.              */
            z = rssringoccs_CDouble_Zero;

            /*  Start the iteration process.                                  */
            for(n = 0; n < maxIterations; n++)
            {

                /*  Calculate real and imaginary parts.                       */
                z = square(z, c);

                /*  Check for divergence.                                     */
                norm = rssringoccs_CDouble_Abs(z);

                if(norm > radius)
                    break;
            }

            if(n == maxIterations)
            {
                red = 0;
                green = 0;
                blue = 0;
            }
            else if (n < 64)
            {
                red = (unsigned char)4*n;
                green = red;
                blue = 255-red;
            }
            else
            {
                red = 255;
                green = 255;
                blue = 0;
            }

            rssringoccs_Color(red, green, blue, fp);
        }
    }
    return 0;
}
