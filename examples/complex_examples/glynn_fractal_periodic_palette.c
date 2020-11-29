#include <stdio.h>
#include <math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

static rssringoccs_ComplexDouble
glynn_func(rssringoccs_ComplexDouble z, rssringoccs_ComplexDouble e,
           rssringoccs_ComplexDouble mu)
{
    rssringoccs_ComplexDouble out;


    out = rssringoccs_Complex_Pow(z, e);
    out = rssringoccs_Complex_Add(out, mu);
    return out;
}

int main(void)
{
    /*  Declare a variable for the output file and give it write permission.  */
    FILE *fp;
    unsigned int x, y, n;
    double z_x, z_y, norm;
    double rcp_factor;
    int maxIterations = 256;
    rssringoccs_ComplexDouble z, e, mu;

    unsigned int size = 1024;
    double red, green, blue;
    unsigned char r, g, b;

    const double x_min = 0.065;
    const double x_max = 0.425;
    const double y_min = -0.67;
    const double y_max = -0.31;

    double radius = 4.0;

    e  = rssringoccs_Complex_Rect(1.5, 0.0);
    mu = rssringoccs_Complex_Rect(-0.2, 0.0);

    fp = fopen("glynn_fractal_periodic_palette.ppm", "w");
    fprintf(fp, "P6\n%d %d\n255\n", size, size);

    rcp_factor = 1.0/(size-1.0);

    /*  Loop through each pixel.                                              */
    for (y=0; y<size; ++y)
    {
        /*  We want to center z_y so scale and shift. This makes the output   *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        z_y = (y - 1.0) * (y_max - y_min) * rcp_factor + y_min;

        for (x=0; x<size; ++x)
        {
            /*  Similarly, center z_x.                                        */
            z_x = x * (x_max - x_min) * rcp_factor + x_min;

            /*  Compute the complex number z_x + i z_y.                       */
            z = rssringoccs_Complex_Rect(z_x, z_y);

            /*  Start the iteration process.                                  */
            for(n = 0; n < maxIterations; ++n)
            {

                /*  Calculate real and imaginary parts.                       */
                z = glynn_func(z, e, mu);

                /*  Check for divergence.                                     */
                norm = rssringoccs_Complex_Abs(z);

                if(norm > radius)
                    break;
            }

            red = sin(0.1*n);
            red = red*red;

            green = sin(0.2*n - 0.78);
            green = green*green;

            blue = sin(0.03*n - 1.78);
            blue = blue*blue;

            r = (unsigned char)(255*red);
            g = (unsigned char)(255*green);
            b = (unsigned char)(255*blue);

            rssringoccs_Color(r, g, b, fp);
        }
        printf("%d\n", y);
    }
    return 0;
}
