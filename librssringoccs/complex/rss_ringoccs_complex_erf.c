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
 ******************************************************************************
 *                        rss_ringoccs_complex_erf                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex error function.                  *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CFloat_Erf:                                               *
 *      rssringoccs_CDouble_Erf:                                              *
 *      rssringoccs_CLDouble_Erf:                                             *
 *  Purpose:                                                                  *
 *      Computes the error function of a complex value z.                     *
 *                                                                            *
 *                               z                                            *
 *                               -                                            *
 *                      2       | |                                           *
 *          Erf(z) = -------    |   exp(-t^2) dt                              *
 *                   sqrt(pi) | |                                             *
 *                             -                                              *
 *                             0                                              *
 *  Arguments:                                                                *
 *      z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          A complex number.                                                 *
 *      w (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          Another complex number.                                           *
 *  Output:                                                                   *
 *      erf_z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):     *
 *          The quotient z / w.                                               *
 *  Method:                                                                   *
 *      Use the complementary error function Erfc(z) and apply the formula:   *
 *                                                                            *
 *          Erf(z) = 1 - Erfc(z)                                              *
 *                                                                            *
 *  NOTES:                                                                    *
 *      No actual float or long double algorithms have been implemented by    *
 *      rss_ringoccs. The complementary error functions simply convert float  *
 *      and long double inputs to doubles
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code is compatible with the C89/C90 standard. The setup script that  *
 *  is used to compile this in config_librssringoccs.sh uses gcc and has the  *
 *  -pedantic and -std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 23, 2020                                             *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Single precision complex error function.                                  */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Erf(rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat erfz, erfcz;

    /*  Compute the complementary error function Erfc(z).                     */
    erfcz = rssringoccs_CFloat_Erfc(z);

    /*  Erf(z) = 1 - Erfc(z). Compute this and return.                        */
    erfz = rssringoccs_CFloat_Subtract_Real(1.0F, erfcz);
    return erfz;
}
/*  End of rssringoccs_CFloat_Erf.                                            */

/*  Double precision complex error function.                                  */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Erf(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble erfz, erfcz;

    /*  Compute the complementary error function Erfc(z).                     */
    erfcz = rssringoccs_CDouble_Erfc(z);

    /*  Erf(z) = 1 - Erfc(z). Compute this and return.                        */
    erfz = rssringoccs_CDouble_Subtract_Real(1.0, erfcz);
    return erfz;
}
/*  End of rssringoccs_CDouble_Erf.                                           */

/*  Double precision complex error function.                                  */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Erf(rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble erfz, erfcz;

    /*  Compute the complementary error function Erfc(z).                     */
    erfcz = rssringoccs_CLDouble_Erfc(z);

    /*  Erf(z) = 1 - Erfc(z). Compute this and return.                        */
    erfz = rssringoccs_CLDouble_Subtract_Real(1.0L, erfcz);
    return erfz;
}
/*  End of rssringoccs_CDouble_Erf.                                           */
