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

/*  The goal of this file is to show off what rss_ringoccs can do and create  *
 *  a really pretty picture. We'll use the complex functions, numerical       *
 *  routines, and the geometry tools provided to wrap the Newton fractal up   *
 *  onto a sphere.                                                            */

/*  We'll need stdio to write to the output ppm file.                         */
#include <stdio.h>

/*  Needed for the "exit" function.                                           */
#include <stdlib.h>

/*  Complex numbers and routines here.                                        */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  2D and 3D vectors are defined here, as well as the stereographic and      *
 *  inverse orthographic projections.                                         */
#include <rss_ringoccs/include/rss_ringoccs_geometry.h>

/*  Routines for creating the PPM figure found here.                          */
#include <rss_ringoccs/include/rss_ringoccs_ppm_plot.h>

static double coeffsa[32] = {
    0.0,  341.0,    0.0, 0.0, 0.0, 0.0,  67518.0,   0.0, 0.0, 0.0,
    0.0, -398505.0, 0.0, 0.0, 0.0, 0.0, -1060200.0, 0.0, 0.0, 0.0,
    0.0,  326895.0, 0.0, 0.0, 0.0, 0.0,  10602.0,   0.0, 0.0, 0.0,
    0.0, -19.0
};

static double coeffsb[31] = {
    -19.0,     0.0, 0.0, 0.0, 0.0, -10602.0,   0.0, 0.0, 0.0, 0.0,
     326895.0, 0.0, 0.0, 0.0, 0.0,  1060200.0, 0.0, 0.0, 0.0, 0.0,
    -398505.0, 0.0, 0.0, 0.0, 0.0, -67518.0,   0.0, 0.0, 0.0, 0.0,
     341.0
};

static rssringoccs_ComplexDouble func(rssringoccs_ComplexDouble z)
{
    rssringoccs_ComplexDouble f, g, out;

    f = rssringoccs_CDouble_Poly_Real_Coeffs(coeffsa, 31, z);
    g = rssringoccs_CDouble_Poly_Real_Coeffs(coeffsb, 30, z);
    out = rssringoccs_CDouble_Divide(f, g);
    return out;
}

/*  This variable acts as the location of the observer looking at the sphere. *
 *  we'll perform our inverse orthographic projection about this point. Note  *
 *  that this variable will be normalize so u and lambda*u will result in the *
 *  same drawing for all positive lambda. u and -u will produce drawings of   *
 *  the opposite side of the sphere.                                          */
static const rssringoccs_ThreeVector camera_pos = {{1.0, 1.0, 1.0}};

/*  The number of pixels in the x and y axes. If you want a higher resolution *
 *  for the output fractal, increase this number. It is best to make n*1024   *
 *  where n is some positive integer.                                         */
static const unsigned int size = 4*1024;

static const unsigned int max_iters = 12;

/*  Tolerance for finding the period of the Julia set.                        */
static const double eps = 1.0e-8;

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
    double z_x, z_y, Pz, norm;
    rssringoccs_ComplexDouble z_0, z_1, z_2, dist;
    double rcp_factor = 1.0 / (size-1.0);

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
    fp = fopen("julia_set_on_sphere.ppm", "w");

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
                    z_0 = *(rssringoccs_ComplexDouble *)(&proj_P);

                    /*  Use the Newton-Raphson function to compute a root.    */
                    for(n = 0; n < max_iters; n++)
                    {
                        z_1 = func(z_0);
                        z_2 = func(z_1);
                        dist = rssringoccs_CDouble_Subtract(z_0, z_2);
                        norm = rssringoccs_CDouble_Abs(dist);

                        if(norm < eps)
                            break;
                        else
                            z_0 = z_1;
                    }
                    if(n == max_iters)
                        rssringoccs_Color(0, 0, 0, fp);
                    else
                        rssringoccs_RGB_Linear_Gradient(n, 0, max_iters, fp);
                }
            }
        }
        printf("%d\n", y);
    }
    return 0;
}
