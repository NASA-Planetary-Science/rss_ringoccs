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
 *                        rss_ringoccs_complex_exp                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Contains the source code for the complex exponention function.        *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_CFloat_Exp:                                               *
 *      rssringoccs_CDouble_Exp:                                              *
 *      rssringoccs_CLDouble_Exp:                                             *
 *  Purpose:                                                                  *
 *      Computes the complex exponential of a complex number.                 *
 *  Arguments:                                                                *
 *      z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):         *
 *          A complex number.                                                 *
 *  Output:                                                                   *
 *      exp_z (rssringoccs_ComplexFloat/ComplexDouble/ComplexLongDouble):     *
 *          The complex exponetial of z.                                      *
 *  Method:                                                                   *
 *      Use Euler's formula. Given z = x + iy we have the following:          *
 *                                                                            *
 *          exp(z) = exp(x + iy) = exp(x) exp(iy)                             *
 *                               = exp*x) [cos(y) + i sin(y)]                 *
 *                               = [exp(x) cos(y)] + i [exp(x) sin(y)]        *
 *                                                                            *
 *      Hence the real and imaginary parts can be computed separately using   *
 *      the real valued exponential and trig functions.                       *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) rss_ringoccs_math.h:                                                  *
 *          This file provides compatibility between the two standard math.h  *
 *          header files (C89 vs C99 math.h). If C99 math.h exists, it simply *
 *          provides aliases for the functions, and if C89 math.h is used     *
 *          it defines the functions missing in the earlier version.          *
 *  2.) rss_ringoccs_complex.h:                                               *
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
 *  Date:       November 12, 2020                                             *
 ******************************************************************************
 *                             Revision History                               *
 ******************************************************************************
 *  2020/12/27 (Ryan Maguire):                                                *
 *      Added float and long double support.                                  *
 *      Frozen for v1.3.                                                      *
 ******************************************************************************/

/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Single precision complex exponential function.                            */
rssringoccs_ComplexFloat
rssringoccs_CFloat_Exp(rssringoccs_ComplexFloat z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexFloat exp_z;
    float real, imag;
    float exp_real, exp_z_real, exp_z_imag;

    /*  Extract the real and imaginary part from z.                           */
    real = rssringoccs_CFloat_Real_Part(z);
    imag = rssringoccs_CFloat_Imag_Part(z);

    /*  We'll use the fact that exp(x+iy) = exp(x)*exp(iy). Then we'll use    *
     *  Euler's formula to write exp(iy) as cos(y) + i*sin(y), giving us      *
     *  exp(z) = exp(x)*cos(y) + i*exp(x)*sin(y).                             */
    exp_real = rssringoccs_Float_Exp(real);

    /*  In the case that z is real, use the real valued exponential. This     *
     *  avoid the result of exp(inf) = inf + i nan. The imaginary part of     *
     *  complex exp(inf) will be exp(inf) * sin(0) = inf * 0 which results in *
     *  nan. This if-then statement avoids this.                              */
    if (imag == 0.0F)
    {
        exp_z_real = exp_real;
        exp_z_imag = 0.0F;
    }

    /*  When we have non-zero imaginary part, resort to Euler's formula.      */
    else
    {
        exp_z_real = exp_real * rssringoccs_Float_Cos(imag);
        exp_z_imag = exp_real * rssringoccs_Float_Sin(imag);
    }

    /*  Use rssringoccs_CFloat_Rect to create the output and return.          */
    exp_z = rssringoccs_CFloat_Rect(exp_z_real, exp_z_imag);
    return exp_z;
}
/*  End of rssringoccs_CFloat_Exp.                                           */

/*  Compute the complex exponential of a complex number z = x + iy.           */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Exp(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble exp_z;
    double real, imag;
    double exp_real, exp_z_real, exp_z_imag;

    /*  Extract the real and imaginary part from z.                           */
    real = rssringoccs_CDouble_Real_Part(z);
    imag = rssringoccs_CDouble_Imag_Part(z);

    /*  We'll use the fact that exp(x+iy) = exp(x)*exp(iy). Then we'll use    *
     *  Euler's formula to write exp(iy) as cos(y) + i*sin(y), giving us      *
     *  exp(z) = exp(x)*cos(y) + i*exp(x)*sin(y).                             */
    exp_real = rssringoccs_Double_Exp(real);

    /*  In the case that z is real, use the real valued exponential. This     *
     *  avoid the result of exp(inf) = inf + i nan. The imaginary part of     *
     *  complex exp(inf) will be exp(inf) * sin(0) = inf * 0 which results in *
     *  nan. This if-then statement avoids this.                              */
    if (imag == 0.0)
    {
        exp_z_real = exp_real;
        exp_z_imag = 0.0;
    }

    /*  When we have non-zero imaginary part, resort to Euler's formula.      */
    else
    {
        exp_z_real = exp_real * rssringoccs_Double_Cos(imag);
        exp_z_imag = exp_real * rssringoccs_Double_Sin(imag);
    }

    /*  Use rssringoccs_CDouble_Rect to create the output and return.         */
    exp_z = rssringoccs_CDouble_Rect(exp_z_real, exp_z_imag);
    return exp_z;
}
/*  End of rssringoccs_CDouble_Exp.                                           */

/*  Long double precision complex exponential.                                */
rssringoccs_ComplexLongDouble
rssringoccs_CLDouble_Exp(rssringoccs_ComplexLongDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexLongDouble exp_z;
    long double real, imag;
    long double exp_real, exp_z_real, exp_z_imag;

    /*  Extract the real and imaginary part from z.                           */
    real = rssringoccs_CLDouble_Real_Part(z);
    imag = rssringoccs_CLDouble_Imag_Part(z);

    /*  We'll use the fact that exp(x+iy) = exp(x)*exp(iy). Then we'll use    *
     *  Euler's formula to write exp(iy) as cos(y) + i*sin(y), giving us      *
     *  exp(z) = exp(x)*cos(y) + i*exp(x)*sin(y).                             */
    exp_real = rssringoccs_LDouble_Exp(real);

    /*  In the case that z is real, use the real valued exponential. This     *
     *  avoid the result of exp(inf) = inf + i nan. The imaginary part of     *
     *  complex exp(inf) will be exp(inf) * sin(0) = inf * 0 which results in *
     *  nan. This if-then statement avoids this.                              */
    if (imag == 0.0L)
    {
        exp_z_real = exp_real;
        exp_z_imag = 0.0L;
    }

    /*  When we have non-zero imaginary part, resort to Euler's formula.      */
    else
    {
        exp_z_real = exp_real * rssringoccs_LDouble_Cos(imag);
        exp_z_imag = exp_real * rssringoccs_LDouble_Sin(imag);
    }

    /*  Use rssringoccs_CDouble_Rect to create the output and return.         */
    exp_z = rssringoccs_CLDouble_Rect(exp_z_real, exp_z_imag);
    return exp_z;
}
/*  End of rssringoccs_CLDouble_Exp.                                          */
