/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble
rssringoccs_Complex_Erf(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble erfz, temp;

    temp = rssringoccs_Complex_Erfc(z);
    erfz = rssringoccs_Complex_Subtract(rssringoccs_Complex_One, temp);

    return erfz;
}
