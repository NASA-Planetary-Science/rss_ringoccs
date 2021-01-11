/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble
rssringoccs_CDouble_Sqrt(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double sqrt_r, theta, real, imag;
    rssringoccs_ComplexDouble sqrt_z;

    /*  We need z in polar coordinates, so compute r and theta.               */
    sqrt_r = rssringoccs_Double_Sqrt(rssringoccs_CDouble_Abs(z));
    theta = rssringoccs_CDouble_Argument(z);

    /*  Once in the form r*exp(i*theta), the square root is compute as        *
     *  sqrt(z) = sqrt(r) * exp(i*theta / 2). r is non-negative, so this is   *
     *  well defined for all z.                                               */
    real = sqrt_r*rssringoccs_Double_Cos(0.5*theta);
    imag = sqrt_r*rssringoccs_Double_Sin(0.5*theta);

    /*  Use rssringoccs_ComplexRect to compute and return sqrt_z.             */
    sqrt_z = rssringoccs_CDouble_Rect(real, imag);
    return sqrt_z;
}
