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
 *              rss_ringoccs_two_slit_fraunhofer_diffraction                  *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the Fraunhofer diffraction modeling of   *
 *      a double slit.                                                        *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
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

/*  The Fresnel integrals are found here.                                     */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

/*  Header file containing the prototypes for the functions.                  */
#include <rss_ringoccs/include/rss_ringoccs_diffraction.h>

float
rssringoccs_Float_Two_Slit_Fraunhofer_Diffraction(float x, float z, float a,
                                                  float d, float lambda)
{
    /*  Declare all necessary variables.                                      */
    float sin_theta, var_1, var_2, var_3, norm, scaled_a, scaled_d, out;

    /*  Compute the length of the point (x, z) in the plane.                  */
    norm = rssringoccs_Float_Sqrt(x*x + z*z);

    /*  If norm is zero, the result is undefined. Return NaN in this case.    */
    if (norm == 0.0)
        out = rssringoccs_NaN;
    else
    {
        /*  Using the fact that sin(theta) = x/r, where r is the hypotenus,   *
         *  we can write sin_theta as follows:                                */
        sin_theta = x/norm;

        /*  Scale the input parameters a and b by the reciprocal of the       *
         *  wavelength lambda.                                                */
        scaled_a = a/lambda;
        scaled_d = d/lambda;

        /*  The various factors are in terms of the squares sinc function and *
         *  the ordinary sine function. Since sinc and sine are more          *
         *  expensive than multiplication (computationally) it is better to   *
         *  compute sinc/sine once, and then square this.                     */
        var_1  = rssringoccs_Float_Sinc(scaled_a*sin_theta);
        var_1 *= var_1;

        var_2  = rssringoccs_Float_Sin(TWO_PI*scaled_d*sin_theta);
        var_2 *= var_2;

        var_3  = 2.0*rssringoccs_Float_Sin(ONE_PI*scaled_d*sin_theta);
        var_3 *= var_3;

        out = var_1*var_2/var_3;
    }

    return out;
}

double
rssringoccs_Double_Two_Slit_Fraunhofer_Diffraction(double x, double z, double a,
                                                   double d, double lambda)
{
    /*  Declare all necessary variables.                                      */
    double sin_theta, var_1, var_2, var_3, norm, scaled_a, scaled_d, out;

    /*  Compute the length of the point (x, z) in the plane.                  */
    norm = rssringoccs_Double_Sqrt(x*x + z*z);

    /*  If norm is zero, the result is undefined. Return NaN in this case.    */
    if (norm == 0.0)
        out = rssringoccs_NaN;
    else
    {
        /*  Using the fact that sin(theta) = x/r, where r is the hypotenus,   *
         *  we can write sin_theta as follows:                                */
        sin_theta = x/norm;

        /*  Scale the input parameters a and b by the reciprocal of the       *
         *  wavelength lambda.                                                */
        scaled_a = a/lambda;
        scaled_d = d/lambda;

        /*  The various factors are in terms of the squares sinc function and *
         *  the ordinary sine function. Since sinc and sine are more          *
         *  expensive than multiplication (computationally) it is better to   *
         *  compute sinc/sine once, and then square this.                     */
        var_1  = rssringoccs_Double_Sinc(scaled_a*sin_theta);
        var_1 *= var_1;

        var_2  = rssringoccs_Double_Sin(TWO_PI*scaled_d*sin_theta);
        var_2 *= var_2;

        var_3  = 2.0*rssringoccs_Double_Sin(ONE_PI*scaled_d*sin_theta);
        var_3 *= var_3;

        out = var_1*var_2/var_3;
    }

    return out;
}

long double
rssringoccs_LDouble_Two_Slit_Fraunhofer_Diffraction(long double x,
                                                       long double z,
                                                       long double a,
                                                       long double d,
                                                       long double lambda)
{
    /*  Declare all necessary variables.                                      */
    long double sin_theta, var_1, var_2, var_3, norm, scaled_a, scaled_d, out;

    /*  Compute the length of the point (x, z) in the plane.                  */
    norm = rssringoccs_LDouble_Sqrt(x*x + z*z);

    /*  If norm is zero, the result is undefined. Return NaN in this case.    */
    if (norm == 0.0)
        out = rssringoccs_NaN;
    else
    {
        /*  Using the fact that sin(theta) = x/r, where r is the hypotenus,   *
         *  we can write sin_theta as follows:                                */
        sin_theta = x/norm;

        /*  Scale the input parameters a and b by the reciprocal of the       *
         *  wavelength lambda.                                                */
        scaled_a = a/lambda;
        scaled_d = d/lambda;

        /*  The various factors are in terms of the squares sinc function and *
         *  the ordinary sine function. Since sinc and sine are more          *
         *  expensive than multiplication (computationally) it is better to   *
         *  compute sinc/sine once, and then square this.                     */
        var_1  = rssringoccs_LDouble_Sinc(scaled_a*sin_theta);
        var_1 *= var_1;

        var_2  = rssringoccs_LDouble_Sin(TWO_PI*scaled_d*sin_theta);
        var_2 *= var_2;

        var_3  = 2.0*rssringoccs_LDouble_Sin(ONE_PI*scaled_d*sin_theta);
        var_3 *= var_3;

        out = var_1*var_2/var_3;
    }

    return out;
}
