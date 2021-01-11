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
 *              rss_ringoccs_one_slit_fraunhofer_diffraction                  *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the Fraunhofer diffraction modeling of   *
 *      a single slit.                                                        *
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
rssringoccs_Float_One_Slit_Fraunhofer_Diffraction(float x, float z, float a)
{
    /*  Declare necessary variables.                                          */
    float result;

    /*  If z is zero we'll return 0.0. Sinc(x) tends to zero for both         *
     *  x -> +infinity and x -> -infinity.                                    */
    if (z == 0.0F)
        result = 0.0F;
    else
    {
        /*  Single slit is computed in terms of the sinc function.            */
        result  = rssringoccs_Float_Sinc(a*x/z);
        result *= result;
    }

    return result;
}

double
rssringoccs_Double_One_Slit_Fraunhofer_Diffraction(double x, double z, double a)
{
    /*  Declare necessary variables.                                          */
    double result;

    /*  If z is zero we'll return 0.0. Sinc(x) tends to zero for both         *
     *  x -> +infinity and x -> -infinity.                                    */
    if (z == 0.0)
        result = 0.0;
    else
    {
        /*  Single slit is computed in terms of the sinc function.            */
        result  = rssringoccs_Double_Sinc(a*x/z);
        result *= result;
    }

    return result;
}

long double
rssringoccs_LDouble_One_Slit_Fraunhofer_Diffraction(long double x,
                                                       long double z,
                                                       long double a)
{
    /*  Declare necessary variables.                                          */
    long double result;

    /*  If z is zero we'll return 0.0. Sinc(x) tends to zero for both         *
     *  x -> +infinity and x -> -infinity.                                    */
    if (z == 0.0)
        result = 0.0;
    else
    {
        /*  Single slit is computed in terms of the sinc function.            */
        result  = rssringoccs_LDouble_Sinc(a*x/z);
        result *= result;
    }

    return result;
}
