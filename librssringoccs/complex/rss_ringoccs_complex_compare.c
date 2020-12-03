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
 *                     rss_ringoccs_complex_compare                           *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for comparing complex data types.            *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 30, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/02 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Boolean data type and True/False are defined here.                        */
#include <rss_ringoccs/include/rss_ringoccs_bool.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  C does not allow structures to be compared, so we need to compare the     *
 *  members of the two complex structs.                                       */

/*  Single precision complex comparison.                                      */
rssringoccs_Bool
rssringoccs_CFloat_Compare(rssringoccs_ComplexFloat z,
                           rssringoccs_ComplexFloat w)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    float z_real, z_imag, w_real, w_imag;

    /*  Extract the real and imaginary parts from z and w.                    */
    z_real = rssringoccs_CFloat_Real_Part(z);
    z_imag = rssringoccs_CFloat_Imag_Part(z);
    w_real = rssringoccs_CFloat_Real_Part(w);
    w_imag = rssringoccs_CFloat_Imag_Part(w);

    /*  Two complex numbers are equal if and only if their real and imaginary *
     *  parts are both equal, so check this.                                  */
    if ((z_real == w_real) && (z_imag == w_imag))
        return rssringoccs_True;
    else
        return rssringoccs_False;
}
/*  End of rssringoccs_CFloat_Compare.                                        */

/*  Double precision complex comparison.                                      */
rssringoccs_Bool
rssringoccs_CDouble_Compare(rssringoccs_ComplexDouble z,
                            rssringoccs_ComplexDouble w)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double z_real, z_imag, w_real, w_imag;

    /*  Extract the real and imaginary parts from z and w.                    */
    z_real = rssringoccs_CDouble_Real_Part(z);
    z_imag = rssringoccs_CDouble_Imag_Part(z);
    w_real = rssringoccs_CDouble_Real_Part(w);
    w_imag = rssringoccs_CDouble_Imag_Part(w);

    /*  Two complex numbers are equal if and only if their real and imaginary *
     *  parts are both equal, so check this.                                  */
    if ((z_real == w_real) && (z_imag == w_imag))
        return rssringoccs_True;
    else
        return rssringoccs_False;
}
/*  End of rssringoccs_CDouble_Compare.                                       */

/*  Long double precision complex comparison.                                 */
rssringoccs_Bool
rssringoccs_CLDouble_Compare(rssringoccs_ComplexLongDouble z,
                             rssringoccs_ComplexLongDouble w)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    long double z_real, z_imag, w_real, w_imag;

    /*  Extract the real and imaginary parts from z and w.                    */
    z_real = rssringoccs_CLDouble_Real_Part(z);
    z_imag = rssringoccs_CLDouble_Imag_Part(z);
    w_real = rssringoccs_CLDouble_Real_Part(w);
    w_imag = rssringoccs_CLDouble_Imag_Part(w);

    /*  Two complex numbers are equal if and only if their real and imaginary *
     *  parts are both equal, so check this.                                  */
    if ((z_real == w_real) && (z_imag == w_imag))
        return rssringoccs_True;
    else
        return rssringoccs_False;
}
/*  End of rssringoccs_CLDouble_Compare.                                      */

#else
/*  Else statement for #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.               */

/*  C99 allows complex values to be compared just like real numbers.          */

/*  Single precision complex comparison.                                      */
rssringoccs_Bool
rssringoccs_CFloat_Compare(rssringoccs_ComplexFloat z,
                           rssringoccs_ComplexFloat w)
{
    if (z == w)
        return rssringoccs_True;
    else
        return rssringoccs_False;
}
/*  End of rssringoccs_CFloat_Compare.                                        */

/*  Double precision complex comparison.                                      */
rssringoccs_Bool
rssringoccs_CDouble_Compare(rssringoccs_ComplexDouble z,
                            rssringoccs_ComplexDouble w)
{
    if (z == w)
        return rssringoccs_True;
    else
        return rssringoccs_False;
}
/*  End of rssringoccs_CDouble_Compare.                                       */

/*  Long double precision complex comparison.                                 */
rssringoccs_Bool
rssringoccs_CLDouble_Compare(rssringoccs_ComplexLongDouble z,
                             rssringoccs_ComplexLongDouble w)
{
    if (z == w)
        return rssringoccs_True;
    else
        return rssringoccs_False;
}
/*  End of rssringoccs_CLDouble_Compare.                                      */

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0.                           */
