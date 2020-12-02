/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble rssringoccs_Complex_Tan(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sin_z, cos_z, tan_z;

    /*  Compute sin(z), cos(z), and then return sin(z)/cos(z).                */
    sin_z = rssringoccs_CDouble_Sin(z);
    cos_z = rssringoccs_CDouble_Cos(z);
    tan_z = rssringoccs_CDouble_Divide(sin_z, cos_z);
    return tan_z;
}
