/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/src/math/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include "rss_ringoccs_complex.h"

rssringoccs_ComplexDouble rssringoccs_Complex_Exp(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble exp_z;
    double real, imag;
    double exp_real, exp_z_real, exp_z_imag;

    /*  Extract the real and imaginary part from z.                           */
    real = rssringoccs_Complex_Real_Part(z);
    imag = rssringoccs_Complex_Imag_Part(z);

    /*  We'll use the fact that exp(x+iy) = exp(x)*exp(iy). Then we'll use    *
     *  Euler's formula to write exp(iy) as cos(y) + i*sin(y), giving us      *
     *  exp(z) = exp(x)*cos(y) + i*exp(x)*sin(y).                             */
    exp_real = rssringoccs_Exp_Double(real);
    exp_z_real = exp_real * rssringoccs_Cos_Double(imag);
    exp_z_imag = exp_real * rssringoccs_Sin_Double(imag);

    /*  Use rssringoccs_ComplexRect to create the output and return.          */
    exp_z = rssringoccs_Complex_Rect(exp_z_real, exp_z_imag);
    return exp_z;
}

rssringoccs_ComplexDouble
rssringoccs_Complex_Log(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double r, theta, real;
    rssringoccs_ComplexDouble ln_z;

    /*  Get the polar representation of the complex number z.                 */
    r = rssringoccs_Complex_Abs(z);
    theta = rssringoccs_Complex_Argument(z);

    /*  The real part is just ln(r), and the imaginary part is theta.         */
    real = rssringoccs_Log_Double(r);

    /*  Use rssringoccs_Complex_Rect to create the complex number and return. */
    ln_z = rssringoccs_Complex_Rect(real, theta);
    return ln_z;
}
