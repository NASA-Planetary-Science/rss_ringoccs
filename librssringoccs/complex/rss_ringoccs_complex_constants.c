/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If C99 complex.h was included before rss_ringoccs_complex.h then we'll    *
 *  use the built-in +, *, -, and / for complex arithmetic.                   */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 1

const rssringoccs_ComplexDouble rssringoccs_Imaginary_Unit = _Complex_I;
const rssringoccs_ComplexDouble rssringoccs_Complex_Zero = 0.0;
const rssringoccs_ComplexDouble rssringoccs_Complex_One = 1.0;
const rssringoccs_ComplexDouble rssringoccs_Complex_NaN = rssringoccs_NaN;
const rssringoccs_ComplexDouble
    rssringoccs_Complex_Infinity = rssringoccs_Infinity +
                                   _Complex_I*rssringoccs_Infinity;

/*  If you do not have complex.h supported, we'll need to define complex      *
 *  arithmetic using the rssringoccs_ComplexDouble data type.                 */
#else

/*  Useful constants.                                                         */
const rssringoccs_ComplexDouble rssringoccs_Imaginary_Unit = {{0.0, 1.0}};
const rssringoccs_ComplexDouble rssringoccs_Complex_Zero   = {{0.0, 0.0}};
const rssringoccs_ComplexDouble rssringoccs_Complex_One    = {{1.0, 0.0}};
const rssringoccs_ComplexDouble
    rssringoccs_Complex_NaN = {{rssringoccs_NaN, rssringoccs_NaN}};
const rssringoccs_ComplexDouble
    rssringoccs_Complex_Infinity = {{rssringoccs_Infinity,
                                     rssringoccs_Infinity}};


#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 1                            */
