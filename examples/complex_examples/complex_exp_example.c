/*  Let's compute the complex exponential of the values 1, i pi, and i.       */

/*  If rss_ringoccs built correctly, rss_ringoccs_complex.h is located in     *
 *  /usr/local/include/rss_ringoccs/include. We can include this as follows:  */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  We'll use stdio to print the results.                                     */
#include <stdio.h>

/*  Define a constant for the value pi.                                       */
const double pi = 3.14159265358979323846264338327950288419716939937510;

int main(void)
{
    /*  Declare all necessary variables. C89 requires declarations at the top *
     *  of a block.                                                           */
    rssringoccs_ComplexDouble z[3];
    rssringoccs_ComplexDouble w[3];
    double re_z, im_z, re_w, im_w;

    /*  And declare a variable for indexing.                                  */
    int n;

    /*  Set z0, z1, and z2 to 1, i, and i pi, respectively.                   */
    z[0] = rssringoccs_Complex_One;
    z[1] = rssringoccs_Imaginary_Unit;
    z[2] = rssringoccs_Complex_Rect(0.0, pi);

    /*  Loop over the results and print them.                                 */
    for (n=0; n<3; ++n)
    {
        /*  Compute the complex exponential of the nth value.                 */
        w[n] = rssringoccs_Complex_Exp(z[n]);

        /*  Extract the real and imaginary parts from z[n] and w[n].          */
        re_z = rssringoccs_Complex_Real_Part(z[n]);
        im_z = rssringoccs_Complex_Imag_Part(z[n]);
        re_w = rssringoccs_Complex_Real_Part(w[n]);
        im_w = rssringoccs_Complex_Imag_Part(w[n]);

        /*  And finally, print the result to the screen.                      */
        printf("exp(%f + i%f) = %f + i%f\n", re_z, im_z, re_w, im_w);
    }

    return 0;
}

/******************************************************************************
 *  We can compile this with:                                                 *
 *                                                                            *
 *      gcc -L/usr/local/lib/ complex_exp_example.c -o test -lrssringoccs     *
 *                                                                            *
 *  If librssringoccs is not in /usr/local/lib/ (this is the default          *
 *  location it is placed in when built via config_src.sh), then change       *
 *  the -L option to the correct location. If /usr/local/include/ is not in   *
 *  your path, add the -I option as follows:                                  *
 *                                                                            *
 *      gcc -I/usr/local/include/ -L/usr/local/lib/                           *
 *              complex_exp_example.c -o test -lrssringoccs                   *
 *                                                                            *
 *  Note, this should all be one line. This outputs an executable "test".     *
 *  Running the executable with ./test, this outputs:                         *
 *      exp(1.000000 + i0.000000) = 2.718282 + i0.000000                      *
 *      exp(0.000000 + i1.000000) = 0.540302 + i0.841471                      *
 *      exp(0.000000 + i3.141593) = -1.000000 + i0.000000                     *
 *  In agreement with known values of the complex exponential.                *
 ******************************************************************************/
