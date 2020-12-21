/*  We'll need stdio to write to the output ppm file.                         */
#include <stdio.h>

/*  Needed for the "exit" function.                                           */
#include <stdlib.h>

/*  Complex numbers and routines here.                                        */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  2D and 3D vectors are defined here, as well as the stereographic and      *
 *  inverse orthographic projections.                                         */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Routines for creating the PPM figure found here.                          */
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

static rssringoccs_ComplexDouble
glynn_func(rssringoccs_ComplexDouble z, rssringoccs_ComplexDouble e,
           rssringoccs_ComplexDouble mu)
{
    rssringoccs_ComplexDouble out;


    out = rssringoccs_CDouble_Pow(z, e);
    out = rssringoccs_CDouble_Add(out, mu);
    return out;
}

/*  This variable acts as the location of the observer looking at the sphere. *
 *  we'll perform our inverse orthographic projection about this point. Note  *
 *  that this variable will be normalize so u and lambda*u will result in the *
 *  same drawing for all positive lambda. u and -u will produce drawings of   *
 *  the opposite side of the sphere.                                          */
static const rssringoccs_ThreeVector camera_pos = {{-0.05, -0.1, -1.0}};

/*  The number of pixels in the x and y axes. If you want a higher resolution *
 *  for the output fractal, increase this number. It is best to make n*1024   *
 *  where n is some positive integer.                                         */
static const unsigned int size = 4*1024;

static const unsigned int max_iters = 256;

/* Values for the min and max of the x and y axes.                            */
static const double x_min = -1.0;
static const double x_max =  1.0;
static const double y_min = -1.0;
static const double y_max =  1.0;

int main(void)
{
    /*  Declare some necessary variables for the geometry.                    */
    rssringoccs_TwoVector proj_P, planar_Z;
    rssringoccs_ThreeVector P, u;

    /*  More dummy variables to loop over.                                    */
    unsigned int x, y, n;
    double z_x, z_y, Pz, norm, red, green, blue;
    unsigned char r, g, b;
    rssringoccs_ComplexDouble z, e, mu;
    double rcp_factor = 1.0 / (size-1.0);

    double radius = 4.0;

    e  = rssringoccs_CDouble_Rect(1.5, 0.0);
    mu = rssringoccs_CDouble_Rect(-0.2, 0.0);

    /*  Declare a variable for the output.                                    */
    FILE *fp;

    /*  Normalize the camera vector and set this to u. First check that the   *
     *  user didn't provide the zero vector since this will cause an error.   */
    norm = rssringoccs_ThreeVector_Euclidean_Norm(camera_pos);

    if (norm == 0.0)
    {
        puts(
            "\nError Encountered: rss_ringoccs\n"
            "\tnewton_fractal_on_sphere.c\n\n"
            "Your input camera_pos is the zero vector. Aborting.\n"
        );
        exit(0);
    }
    else
        u = rssringoccs_ThreeVector_Normalize(camera_pos);

    /* Open and name the file and give it write permission.                   */
    fp = fopen("glynn_fractal_on_sphere.ppm", "w");

    /*  Needed to create the output ppm file. This is the preamble.           */
    fprintf(fp, "P6\n%d %d\n255\n", size, size);

    /*  Loop over every pixel and perform the Newton-Raphson method with the  *
     *  given point as the initial guess for the algorithm.                   */
    for (y=0; y<size; ++y)
    {
        /*  We want to center z_y so scale and shift. This makes the output   *
         *  picture lie in the box [x_min, x_max] x [y_min, y_max].           */
        z_y = (size - y - 1.0) * (y_max - y_min) * rcp_factor + y_min;

        for (x=0; x<size; ++x)
        {
            /*  Similarly, center z_x.                                        */
            z_x = x * (x_max - x_min) * rcp_factor + x_min;

            if ((z_x*z_x + z_y*z_y >= 1.0))
                rssringoccs_Color(0, 0, 0, fp);
            else
            {
                planar_Z = rssringoccs_TwoVector_Rect(z_x, z_y);
                P = rssringoccs_Inverse_Orthographic_Projection(planar_Z, u);
                Pz = rssringoccs_ThreeVector_Z(P);

                if (Pz > 0.999999)
                    rssringoccs_Color(128, 128, 128, fp);
                else
                {
                    proj_P = rssringoccs_Stereographic_Projection(P);

                    /*  rssringoccs_ComplexDouble and rssringoccs_TwoVector   *
                     *  define identical structures, so we can just cast      *
                     *  proj_P to a rssringoccs_ComplexDouble pointer and     *
                     *  perform Newton-Raphson with that.                     */
                    z = *(rssringoccs_ComplexDouble *)(&proj_P);

                    /*  Use the Newton-Raphson function to compute a root.    */
                    for(n = 0; n < max_iters; ++n)
                    {

                        /*  Calculate real and imaginary parts.               */
                        z = glynn_func(z, e, mu);

                        /*  Check for divergence.                             */
                        norm = rssringoccs_CDouble_Abs(z);

                        if(norm > radius)
                            break;
                    }
                    red = rssringoccs_Double_Sin(0.1*n);
                    red = red*red;

                    green = rssringoccs_Double_Sin(0.2*n - 0.78);
                    green = green*green;

                    blue = rssringoccs_Double_Sin(0.03*n - 1.78);
                    blue = blue*blue;

                    r = (unsigned char)(255*red);
                    g = (unsigned char)(255*green);
                    b = (unsigned char)(255*blue);

                    rssringoccs_Color(r, g, b, fp);
                }
            }
        }
    }
    return 0;
}
