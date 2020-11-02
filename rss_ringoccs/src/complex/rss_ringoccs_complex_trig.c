/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/src/math/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include "rss_ringoccs_complex.h"

rssringoccs_ComplexDouble rssringoccs_Complex_Sin(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double x, y, real, imag;
    rssringoccs_ComplexDouble sin_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_Complex_Real_Part(z);
    y = rssringoccs_Complex_Imag_Part(z);

    /*  The real part is sin(x)cosh(y).                                       */
    real = rssringoccs_Sin_Double(x)*rssringoccs_Cosh_Double(y);
    imag = rssringoccs_Cos_Double(x)*rssringoccs_Sinh_Double(y);

    /*  Use rssringoccs_Complex_Rect to create the output and return.         */
    sin_z = rssringoccs_Complex_Rect(real, imag);
    return sin_z;
}

rssringoccs_ComplexDouble rssringoccs_Complex_Cos(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double x, y, real, imag;
    rssringoccs_ComplexDouble cos_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_Complex_Real_Part(z);
    y = rssringoccs_Complex_Imag_Part(z);

    /*  The real part is sin(x)cosh(y).                                       */
    real = rssringoccs_Cos_Double(x)*rssringoccs_Cosh_Double(y);
    imag = -rssringoccs_Sin_Double(x)*rssringoccs_Sinh_Double(y);

    /*  Use rssringoccs_Complex_Rect to create the output and return.         */
    cos_z = rssringoccs_Complex_Rect(real, imag);
    return cos_z;
}

rssringoccs_ComplexDouble rssringoccs_Complex_Tan(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sin_z, cos_z, tan_z;

    /*  Compute sin(z), cos(z), and then return sin(z)/cos(z).                */
    sin_z = rssringoccs_Complex_Sin(z);
    cos_z = rssringoccs_Complex_Cos(z);
    tan_z = rssringoccs_Complex_Divide(sin_z, cos_z);
    return tan_z;
}
