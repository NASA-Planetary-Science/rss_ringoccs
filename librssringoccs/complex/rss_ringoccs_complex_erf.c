/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble
rssringoccs_CDouble_Erf(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble erfz, temp;

    temp = rssringoccs_CDouble_Erfc(z);
    erfz = rssringoccs_CDouble_Subtract_Real(1.0, temp);

    return erfz;
}
