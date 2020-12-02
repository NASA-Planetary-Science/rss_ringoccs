/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble
rssringoccs_CDouble_Sin(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double x, y, real, imag;
    rssringoccs_ComplexDouble sin_z;

    /*  Extract the real and imaginary parts from z.                          */
    x = rssringoccs_CDouble_Real_Part(z);
    y = rssringoccs_CDouble_Imag_Part(z);

    /*  The real part is sin(x)cosh(y).                                       */
    real = rssringoccs_Double_Sin(x)*rssringoccs_Double_Cosh(y);
    imag = rssringoccs_Double_Cos(x)*rssringoccs_Double_Sinh(y);

    /*  Use rssringoccs_Complex_Rect to create the output and return.         */
    sin_z = rssringoccs_CDouble_Rect(real, imag);
    return sin_z;
}
