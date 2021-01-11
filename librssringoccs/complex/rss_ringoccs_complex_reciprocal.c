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
 *                        rss_ringoccs_complex_add                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for complex addition.                        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_complex.h:                                               *
 *          Header where complex types and function prototypes are defined.   *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 30, 2020                                             *
 ******************************************************************************/

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If _RSS_RINGOCCS_USING_COMPLEX_H_ is set to zero, then C99 complex.h has  *
 *  not been included and we must define our own algorithms.                  */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 0

/*  In C99, since _Complex is a built-in data type, given double _Complex z1  *
 *  and double _Complex z2, you can just do z1 * z2. Structs cannot be        *
 *  multiplied so we need a function for computing this.                      */

/*  Function for computing the reciprocal (or inverse) of a complex number.   */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Reciprocal(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble conj_z, rcpr_z;
    double rcp_abs_z_sq;

    /*  Compute the conjugate of z and its absolute value.                    */
    conj_z = rssringoccs_CDouble_Conjugate(z);
    rcp_abs_z_sq = rssringoccs_CDouble_Abs_Squared(z);

    /*  The inverse of z is conj_z / abs_z^2, so return this.                 */
    rcpr_z = rssringoccs_CDouble_Multiply_Real(rcp_abs_z_sq, conj_z);
    return rcpr_z;
}

#else

/*  Function for computing 1/z for non-zero z.                                */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Reciprocal(rssringoccs_ComplexDouble z)
{
    return 1.0/z;
}

#endif
