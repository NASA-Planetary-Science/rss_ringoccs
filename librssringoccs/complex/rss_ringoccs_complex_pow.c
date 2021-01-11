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
 *                        rss_ringoccs_complex_pow                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the complex power function. If           *
 *      rss_ringoccs was built with C99 complex.h support, this function is   *
 *      not compiled and instead rssringoccs_Complex_Pow is just an alias for *
 *      the cpow function. By default rss_ringoccs builds with C89 (commonly  *
 *      called ANSI C) support, so this file will be a part of the build.     *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_complex.h:                                               *
 *          Header file where rssringoccs_ComplexDouble is defined, as well   *
 *          as the prototype for rssringoccs_Complex_Pow.                     *
 ******************************************************************************
 *                                 WARNINGS                                   *
 *  1.) This function implicitly uses the complex log function, and hence     *
 *      there is a branch cut in the second variable. The function is         *
 *      continuous in the first variable, i.e. no branch cut.                 *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 12, 2020                                             *
 ******************************************************************************
 *                                History                                     *
 *  2020/11/14 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Compute the complex power z0^z1.                                          */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Pow(rssringoccs_ComplexDouble z0,
                        rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z0_to_the_z1, ln_z0, z1_ln_z0;

    /*  We can write x^y as exp(y ln(x)) and this is how we'll compute for    *
     *  complex powers. First compute log(z1).                                */
    ln_z0 = rssringoccs_CDouble_Log(z0);

    /*  Next use rssringoccs_Complex_Multiply to compute the product with z0. */
    z1_ln_z0 = rssringoccs_CDouble_Multiply(z1, ln_z0);

    /*  And finally exponentiate.                                             */
    z0_to_the_z1 = rssringoccs_CDouble_Exp(z1_ln_z0);
    return z0_to_the_z1;
}

/*  Compute the complex power z^x for x real.                                 */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Real_Pow(rssringoccs_ComplexDouble z, double x)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z_to_the_x, ln_z, x_ln_z;

    /*  We can write z^x as exp(x ln(z)) and this is how we'll compute for    *
     *  complex powers. First compute log(z).                                 */
    ln_z = rssringoccs_CDouble_Log(z);

    /*  Next use rssringoccs_Complex_Scale to compute the product with x.     */
    x_ln_z = rssringoccs_CDouble_Multiply_Real(x, ln_z);

    /*  And finally exponentiate.                                             */
    z_to_the_x = rssringoccs_CDouble_Exp(x_ln_z);
    return z_to_the_x;
}
