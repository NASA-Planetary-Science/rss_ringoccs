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
 *                  rss_ringoccs_gap_fresnel_diffraction                      *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the Fresnel diffraction of a gap.        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_complex.h:                                               *
 *          Header file where rssringoccs_ComplexDouble is defined, as well   *
 *          as the prototype for rssringoccs_Complex_Cos.                     *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       November 27, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/11/27 (Ryan Maguire):                                                *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Definition of rssringoccs_ComplexDouble found here.                       */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  The Fresnel integrals are found here.                                     */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Header file containing the prototypes for the functions.                  */
#include <rss_ringoccs/include/rss_ringoccs_diffraction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      Inverted_Square_Well_Diffraction_Float                                *
 *  Purpose:                                                                  *
 *      Compute the diffraction pattern from a plane wave incident on a       *
 *      square well, assuming the Fresnel approximation is valid.             *
 *  Arguments:                                                                *
 *      x (float):                                                            *
 *          The location on the x-axis for the point being computed.          *
 *      a (float):                                                            *
 *          The left-most endpoint of the square well.                        *
 *      b (float):                                                            *
 *          The right-most endpoint of the square well.                       *
 *      F (float):                                                            *
 *          The Fresnel scale.                                                *
 *  Notes:                                                                    *
 *      1.) This function relies on the C99 standard, or higher.              *
 ******************************************************************************/

rssringoccs_ComplexDouble
rssringoccs_Complex_Gap_Diffraction(double x, double a, double b, double F)
{
    double arg1, arg2;
    rssringoccs_ComplexDouble z1, z2, out, scale;

    /*  The scale factor for the integral is (1-i)/sqrt(2 pi).                */
    scale = rssringoccs_CDouble_Rect(SQRT_ONE_BY_TWO_PI, -SQRT_ONE_BY_TWO_PI);

    /*  The bounds of the integral are sqrt(pi/2)(a-x)/F and                  *
     *  sqrt(pi/2)(b-x)/F, and the output is computed in terms of this.       */
    arg1 = SQRT_PI_BY_2*(a-x)/F;
    arg2 = SQRT_PI_BY_2*(b-x)/F;

    /*  Compute the Fresnel integrals of the two arguments.                   */
    z1 = rssringoccs_Complex_Fresnel_Integral(arg1);
    z2 = rssringoccs_Complex_Fresnel_Integral(arg2);

    /*  The output is the difference of z2 and z1 scaled by the factor        *
     *  sqrt(1/2 pi) * (1 - i), so compute this.                              */
    out = rssringoccs_CDouble_Subtract(z2, z1);
    out = rssringoccs_CDouble_Multiply(out, scale);

    return out;
}
