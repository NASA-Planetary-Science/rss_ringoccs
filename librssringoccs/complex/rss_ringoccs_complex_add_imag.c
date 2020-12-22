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
 *                     rss_ringoccs_complex_add_imag                          *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex addition.                        *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CFloat_Add_Imag:                                          *
 *      rssringoccs_CDouble_Add_Imag:                                         *
 *      rssringoccs_CLDouble_Add_Imag:                                        *
 *  Purpose:                                                                  *
 *      Adds an imaginary number to a complex one.                            *
 *                                                                            *
 *          w = z + iy = (a + ib) + iy = a + i(b+y)                           *
 *                                                                            *
 *  Arguments:                                                                *
 *      y (float/double/long double):                                         *
 *          The imaginary number we wish to add to z.                         *
 *      z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      w (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          The sum of z and iy.                                              *
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
 *  -pedantic and =std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       December 3, 2020                                              *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/03 (Ryan Maguire):                                                *
 *      Moved here from rss_ringoccs_complex_add.c.                           *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  In C99, since _Complex is a built-in data type, doubles and _Complex      *
 *  doubles can be added via y*_Complex_I + z. With C89 we use structs to     *
 *  define complex numbers. Since we can't add a double to a struct, and      *
 *  since the _Complex_I macro is undefined, we need a function for computing *
 *  the sum of complex numbers with imaginary ones.                           */

/*  Single precision complex addition where one variable is imaginary.        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Add_Imag(float y, rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat sum;
    float real, imag, sum_im;

    /*  Extract the real and imaginary parts from the input.                  */
    real = rssringoccs_CFloat_Real_Part(z);
    imag = rssringoccs_CFloat_Imag_Part(z);

    /*  Add the imaginary variable to the imaginary part of z.                */
    sum_im = y + imag;

    /*  Create the output from real and sum_im and return.                    */
    sum = rssringoccs_CFloat_Rect(real, sum_im);
    return sum;
}
/*  End of rssringoccs_CFloat_Add_Imag.                                       */

/*  Double precision complex addition where one variable is imaginary.        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Add_Imag(double y, rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sum;
    double real, imag, sum_im;

    /*  Extract the real and imaginary parts from the input.                  */
    real = rssringoccs_CDouble_Real_Part(z);
    imag = rssringoccs_CDouble_Imag_Part(z);

    /*  Add the imaginary variable to the imaginary part of z.                */
    sum_im = y + imag;

    /*  Create the output from real and sum_im and return.                    */
    sum = rssringoccs_CDouble_Rect(real, sum_im);
    return sum;
}
/*  End of rssringoccs_CDouble_Add_Imag.                                      */

/*  Long double precision complex addition where one variable is imaginary.   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add_Imag(long double y, rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble sum;
    long double real, imag, sum_im;

    /*  Extract the real and imaginary parts from the input.                  */
    real = rssringoccs_CLDouble_Real_Part(z);
    imag = rssringoccs_CLDouble_Imag_Part(z);

    /*  Add the imaginary variable to the imaginary part of z.                */
    sum_im = y + imag;

    /*  Create the output from real and sum_im and return.                    */
    sum = rssringoccs_CLDouble_Rect(real, sum_im);
    return sum;
}
/*  End of rssringoccs_CLDouble_Add_Imag.                                     */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  If we get here we have complex.h support so we'll just use the            *
 *  arithmetic made available in the C99 standard.                            */

/*  Single precision complex addition.                                        */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Add_Imag(float y, rssringoccs_ComplexFloat z)
{
    return (float _Complex)_Complex_I*y + z;
}
/*  End of rssringoccs_CFloat_Add_Imag.                                       */

/*  Double precision complex addition.                                        */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Add_Imag(double y, rssringoccs_ComplexDouble z)
{
    return (double _Complex)_Complex_I*y + z;
}
/*  End of rssringoccs_CDouble_Add_Imag.                                      */

/*  Long double precision complex addition.                                   */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Add_Imag(long double y, rssringoccs_ComplexLongDouble z)
{
    return (long double _Complex)_Complex_I*y + z;
}
/*  End of rssringoccs_CLDouble_Add_Imag.                                     */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */
