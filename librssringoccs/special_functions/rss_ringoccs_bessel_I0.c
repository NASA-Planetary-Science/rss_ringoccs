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
 *                                Bessel I0                                   *
 ******************************************************************************
 *  This file contains functions for computing the Bessel I0 function.        *
 *  We define I_0(x) as follows:                                              *
 *                                                                            *
 *      I_0(x)  =  J_0(ix)                                                    *
 *                                                                            *
 *  Where J_0 is the Bessel function of the First kind with nu = 0 and i is   *
 *  the imaginary unit.                                                       *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  BesselI0:                                                                 *
 *  Purpose:                                                                  *
 *      Compute the I_0 bessel function for a real argument.                  *
 *  Arguments:                                                                *
 *      x (float/double)                                                      *
 *          A real number, the argument for I_0(x).                           *
 *  Output:                                                                   *
 *      bessel_I0:                                                            *
 *          The Bessel function I_0(x).                                       *
 *  Method:                                                                   *
 *      This uses the standard Taylor series for small inputs, and            *
 *      asymptotic expansions for large inputs. The Taylor series is:         *
 *                                                                            *
 *                   ___                                                      *
 *                   \           x^(2n)                                       *
 *        I0(x)=     /__     --------------                                   *
 *                   n = 0   ((n)! * 2^n)^2                                   *
 *                                                                            *
 *      This can be obtained by simply evaluating ix into the J_0 function.   *
 *                                                                            *
 *      The asymptotic expansions is of the form:                             *
 *                                                                            *
 *          I_0(x)  ~   exp(x) / sqrt(2 * pi * x)                             *
 *                                                                            *
 *      For very large values, we simply return infinity.                     *
 *  Error:                                                                    *
 *      In the region in which the Taylor series is used, relative error is   *
 *      10^-16. In the hand-off region with the asymptotic expansion, the     *
 *      error is 10^-9 but quickly drops back to 10^-16.                      *
 *      The regions where the Taylor series is used are listed below:         *
 *          float:           (-3.46, 3.46)                                    *
 *          double:          (-4.00, 4.00)                                    *
 *          long double:     (-4.35, 4.35)                                    *
 *  Notes:                                                                    *
 *      This code was written with an emphasis on clarity and accessibility   *
 *      to a wide audience without the need for anything too sophisticated.   *
 *      More accurate methods involving rational functions and Chebyshev      *
 *      polynomials exist. See the GNU Scientific Library for source code.    *
 ******************************************************************************
 *                            A FRIENDLY WARNING                              *
 ******************************************************************************
 *  This code is compatible with the C89/C90 standard. The setup script that  *
 *  is used to compile this in configspecialfunctions.sh uses gcc and has the *
 *  -pedantic and =std=c89 flags to check for compliance. If you edit this to *
 *  use C99 features (built-in complex, built-in booleans, C++ style comments *
 *  and etc.), or GCC extensions, you will need to edit the config script.    *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) math.h:                                                               *
 *      Standard library for mathematical functions like sin, cos, atan.      *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                          A NOTE ON CONVENTIONS                             *
 ******************************************************************************
 *  1.) i is a complex number, we'll use n for indexing. Max index values for *
 *      a for loop or while loop should be either capital N, or preferrably   *
 *      a variable name describing what this limit is. i.e. if it's the size  *
 *      of an array, call it arr_size.                                        *
 *  2.) Do not cast malloc. While this is required in C++, it is not in C and *
 *      the official documentation even frowns upon it. That is, do this:     *
 *          double *ptr;                                                      *
 *          long N;                                                           *
 *          ...                                                               *
 *          ptr = malloc(sizeof(*ptr) * N);                                   *
 *      and not:                                                              *
 *          ptr = (double *)malloc(sizeof(*ptr) * N);                         *
 *      malloc returns a void pointer and this is automatically and safely    *
 *      promoted to whatever the type of ptr is.                              *
 *  3.) While not frowned upon, the following makes maintenance easier.       *
 *      Instead of using malloc like (for a pointer to a double):             *
 *          ptr = malloc(sizeof(double) * N);                                 *
 *      Use:                                                                  *
 *          ptr = malloc(sizeof(*ptr) * N);                                   *
 *  4.) Declare variables towards the top of the code and not inside a for    *
 *      loop. It's extremely ugly. The following is good:                     *
 *          long k;                                                           *
 *          long N = 100;                                                     *
 *          for (k=0; k<N; ++k)                                               *
 *              do stuff                                                      *
 *      And this looks horrid:                                                *
 *          long N = 100;                                                     *
 *          for (long k=0; k<N; ++k)                                          *
 *              do stuff                                                      *
 *      While legal in C99, it is not in C89 and the way the config script is *
 *      is set up will cause compiler errors.                                 *
 *  5.) If a number is needed in a for loop, even once, declare or set it as  *
 *      a macro, rather than having the number verbatim in the loop. If the   *
 *      the number needs to change it's easier to keep track this way.        *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       September 15, 2019                                            *
 ******************************************************************************
 *                                History                                     *
 ******************************************************************************
 *  2020/09/12 (Ryan Maguire):                                                *
 *      Made code compliant with the C89 standard. Added comments.            *
 ******************************************************************************/

/*  The C standard library header math.h is included here, as are aliases for *
 *  various functions, the macros INFINITY and NAN, as well as the max legal  *
 *  values for the exponential function which don't return INFINITY.          */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Complex variables and functions defined here.                             */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  Prototypes for these functions declared here.                             */
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>

static unsigned int __taylor_float_deg = 16;
static unsigned int __asym_float_deg = 4;
static float __taylor_float[17] = {
    1.0F,                                 0.25F,
    1.56250e-2F,                          4.34027777777777777777777777778e-4F,
    6.78168402777777777777777777778e-6F,  6.78168402777777777777777777778e-8F,
    4.70950279706790123456790123457e-10F, 2.40280754952443940539178634417e-12F,
    9.38596699032984142731166540690e-15F, 2.89690339207711155163940290337e-17F,
    7.24225848019277887909850725841e-20F, 1.49633439673404522295423703686e-22F,
    2.59780277210771740096221707789e-25F, 3.84290350903508491266600159451e-28F,
    4.90166263907536340901275713585e-31F, 5.44629182119484823223639681761e-34F,
    5.31864435663559397679335626720e-37F
};

static float __asym_float[5] = {
    1.0F, 0.1250F, 0.07031250F, 0.07324218750F, 0.1121520996093750F
};

static unsigned int __taylor_double_deg = 24;
static unsigned int __asym_double_deg = 6;
static double __taylor_double[25] = {
    1.0,                                 0.25,
    1.56250e-2,                          4.34027777777777777777777777778e-4,
    6.78168402777777777777777777778e-6,  6.78168402777777777777777777778e-8,
    4.70950279706790123456790123457e-10, 2.40280754952443940539178634417e-12,
    9.38596699032984142731166540690e-15, 2.89690339207711155163940290337e-17,
    7.24225848019277887909850725841e-20, 1.49633439673404522295423703686e-22,
    2.59780277210771740096221707789e-25, 3.84290350903508491266600159451e-28,
    4.90166263907536340901275713585e-31, 5.44629182119484823223639681761e-34,
    5.31864435663559397679335626720e-37, 4.60090342269515049895619054256e-40,
    3.55007980146230748376249270259e-43, 2.45850401763317692781336059736e-46,
    1.53656501102073557988335037335e-49, 8.71068600351890918301219032512e-53,
    4.49932128280935391684513963074e-56, 2.12633330945621640682662553438e-59,
    9.22887721118149482129611777074e-63
};

static double __asym_double[7] = {
    1.0,                        0.1250,
    0.07031250,                 0.07324218750,
    0.1121520996093750,         0.2271080017089843750,
    0.57250142097473144531250
};

static unsigned int __taylor_LDouble_deg = 28;
static unsigned int __asym_LDouble_deg = 8;
static long double __taylor_longdouble[29] = {
    1.0L,                                 0.25L,
    1.56250e-2L,                          4.34027777777777777777777777778e-4L,
    6.78168402777777777777777777778e-6L,  6.78168402777777777777777777778e-8L,
    4.70950279706790123456790123457e-10L, 2.40280754952443940539178634417e-12L,
    9.38596699032984142731166540690e-15L, 2.89690339207711155163940290337e-17L,
    7.24225848019277887909850725841e-20L, 1.49633439673404522295423703686e-22L,
    2.59780277210771740096221707789e-25L, 3.84290350903508491266600159451e-28L,
    4.90166263907536340901275713585e-31L, 5.44629182119484823223639681761e-34L,
    5.31864435663559397679335626720e-37L, 4.60090342269515049895619054256e-40L,
    3.55007980146230748376249270259e-43L, 2.45850401763317692781336059736e-46L,
    1.53656501102073557988335037335e-49L, 8.71068600351890918301219032512e-53L,
    4.49932128280935391684513963074e-56L, 2.12633330945621640682662553438e-59L,
    9.22887721118149482129611777074e-63L, 3.69155088447259792851844710830e-66L,
    1.36521852236412645285445529153e-69L, 4.68181934967121554476836519729e-73L,
    1.49292708854311720177562665730e-76L
};

static long double __asym_longdouble[9] = {
    1.0L,                            0.1250L,
    0.07031250L,                     0.07324218750L,
    0.1121520996093750L,             0.2271080017089843750L,
    0.57250142097473144531250L,      1.72772750258445739746093750L,
    6.07404200127348303794860839844L
};

/*  Compute the Bessel I_0 function for a floating point number x. This       *
 *  returns floating point precision, maximum relative error ~1.e-6.          */
float rssringoccs_Float_Bessel_I0(float x)
{
    /*  Declare necessary variables. C89 requires declaring these at the top. */
    float bessel_I0, arg, abs_x;

    /*  I_0 is an even function, so compute the absolute value of x.          *
     *  rssringoccs_Abs_Float is a macro for fabsf, if your system provides   *
     *  it (C99 math.h) and fabs if not, defined in rss_ringoccs_math.h.      */
    abs_x = rssringoccs_Float_Abs(x);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (abs_x < 12.0F)
    {
        /*  The series is in powers of x^2, so use Horner's method of         *
         *  polynomial computation with that.                                 */
        arg = abs_x*abs_x;

        /*  Compute the Taylor series for bessel I0.                          */
        bessel_I0 = rssringoccs_Real_Poly_Float_Coeffs(
            __taylor_float, __taylor_float_deg, arg
        );
    }

    /*  For larger values, use the asymptotic expansion. MAX_FLOAT_BASE_E is  *
     *  the largest float x which doesn't return INFINITY. It is defined in   *
     *  rss_ringocc_math.h.                                                   */
    else if (abs_x < rssringoccs_Max_Float_Base_E)
    {
        /*  The asymptotic expansion is in terms of 1/x.                      */
        arg = 1.0F/abs_x;

        /*  Compute the polynomial term using Horner's Method.                */
        bessel_I0 = rssringoccs_Real_Poly_Float_Coeffs(
            __asym_float, __asym_float_deg, arg
        );

        /*  Multiply by the coefficient factor and return. Exp_Float and      *
         *  Sqrt_Float are aliases for expf and sqrtf, respectively, if your  *
         *  platform provides it, and exp and sqrt, respectively, otherwise.  *
         *  They are defined in rss_ringoccs_math.h. The rssringoccs_Two_Pi macros is     *
         *  defined in rss_ringoccs_math.h.                                   */
        bessel_I0 *= rssringoccs_Float_Exp(abs_x) /
                     rssringoccs_Float_Sqrt(rssringoccs_Two_Pi_F*abs_x);
    }

    /*  For very large inputs, return INFINITY. INFINITY is standard in C99,  *
     *  but defined in rss_ringoccs_math.h if not available.                  */
    else
        bessel_I0 = rssringoccs_Infinity_F;

    return bessel_I0;
}

/*  Compute the Bessel I_0 function for a double precision value x. This      *
 *  returns double precision, maximum relative error ~1.e-9.                  */
double rssringoccs_Double_Bessel_I0(double x)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    double bessel_I0, arg, abs_x;

    /*  I_0 is symmetric so compute the absolute value of x and use that.     *
     *  rssringoccs_Abs_Double is an alias for fabs.                          */
    abs_x = rssringoccs_Double_Abs(x);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (abs_x < 16.0)
    {
        /*  The series is in powers of x^2, so use Horner's method with that. */
        arg = abs_x*abs_x;

        /*  Compute the Taylor series for bessel I0.                          */
        bessel_I0 = rssringoccs_Real_Poly_Double_Coeffs(
            __taylor_double, __taylor_double_deg, arg
        );
    }

    /*  For larger values, use the asymptotic expansion. MAX_DOUBLE_BASE_E is *
     *  defined in rss_ringoccs_math.h, it's the largest value for a double   *
     *  which won't return INFINITY for exp(x).                               */
    else if (abs_x < rssringoccs_Max_Double_Base_E)
    {
        /*  The asymptotic expansion is in terms of 1/x.                      */
        arg = 1.0/abs_x;

        /*  Compute the polynomial part of the asymptotic expansion.          */
        bessel_I0 = rssringoccs_Real_Poly_Double_Coeffs(
            __asym_double, __asym_double_deg, arg
        );

        /*  Multiply by the coefficient factor and return. Exp_Double and     *
         *  Sqrt_Double are aliases for exp and sqrt, respectively, defined   *
         *  in rss_ringoccs_math.h. rssringoccs_Two_Pi is in rss_ringoccs_math.h.         */
        bessel_I0 *= rssringoccs_Double_Exp(abs_x) /
                     rssringoccs_Double_Sqrt(rssringoccs_Two_Pi*abs_x);
    }

    /*  For very large inputs, return INFINITY. INFINITY is standard in C99,  *
     *  and provided in rss_ringoccs_math.h otherwise.                        */
    else
        bessel_I0 = rssringoccs_Infinity;

    return bessel_I0;
}

/*  Compute the Bessel I_0 function for a long double precision value x. This *
 *  returns long double precision, maximum relative error ~1.e-14.            */
long double rssringoccs_LDouble_Bessel_I0(long double x)
{
    /*  Declare necessary variables.                                          */
    long double bessel_I0, arg, abs_x;

    /*  I_0 is symmetric so compute the absolute value of x and use that.     *
     *  Abs_LongDouble is defined in rss_ringoccs_math.h, it is an alias for  *
     *  fabsl, if available, and fabs otherwise.                              */
    abs_x = rssringoccs_LDouble_Abs(x);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (abs_x < 19.0L)
    {
        /*  The series is in powers of x^2, so use Horner's method with that. */
        arg = abs_x*abs_x;

        bessel_I0 = rssringoccs_Real_Poly_LDouble_Coeffs(
            __taylor_longdouble, __taylor_LDouble_deg, arg
        );
    }

    /*  For larger values, use the asymptotic expansion. MAX_LDOUBLE_BASE_E   *
     *  is defined in rss_ringoccs_math.h, it's the largest value for a       *
     *  long double which won't return INFINITY for exp(x).                   */
    else if (abs_x < rssringoccs_Max_LDouble_Base_E)
    {
        /*  The asymptotic expansion is in terms of 1/x.                      */
        arg = 1.0L/abs_x;

        bessel_I0 = rssringoccs_Real_Poly_LDouble_Coeffs(
            __asym_longdouble, __asym_LDouble_deg, arg
        );

        /*  Multiply by the coefficient factor and return. Exp_LongDouble and *
         *  Sqrt_LongDouble are aliases for expl and sqrtl, respectively, if  *
         *  available, and exp and sqrt otherwise. They are defined in        *
         *  rss_ringoccs_math.h. rssringoccs_Two_Pi is in rss_ringoccs_math.h.            */
        bessel_I0 *= rssringoccs_LDouble_Exp(abs_x) /
                     rssringoccs_LDouble_Sqrt(rssringoccs_Two_Pi_L*abs_x);
    }

    /*  For very large inputs, return INFINITY. This is standard in C99, and  *
     *  provided in rss_ringoccs_math.h otherwise.                            */
    else
        bessel_I0 = rssringoccs_Infinity_L;

    return bessel_I0;
}

/*  Compute the Bessel I_0 function for a complex value z. This               *
 *  returns double precision, maximum error ~1.e-9.                           */
rssringoccs_ComplexDouble
rssringoccs_CDouble_Bessel_I0(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    rssringoccs_ComplexDouble bessel_I0, arg, exp_arg, sqrt_arg;
    double abs_z, abs_z_real, z_real, z_imag, real_bessel;

    /*  Extract the real part of z.                                           */
    z_real = rssringoccs_CDouble_Real_Part(z);

    /*  Compute the magnitude of z and the absolute value of its real part.   */
    abs_z = rssringoccs_CDouble_Abs(z);

    /*  If the real part is zero we obtain the Bessel J0 function.            */
    if (z_real == 0)
    {
        /*  Compute the imaginary part of z and use this to compute Bessel J0.*/
        z_imag = rssringoccs_CDouble_Imag_Part(z);
        real_bessel = rssringoccs_Double_Bessel_J0(z_imag);

        /*  The output is I0 = J0 + 0*Imaginary_Unit.                         */
        bessel_I0 = rssringoccs_CDouble_Rect(real_bessel, 0.0);
        return bessel_I0;
    }

    /*  If z is negative, flip z. That is, negative z so that it lies in the  *
     *  right-half of the complex plane. I0 is even so negation does not      *
     *  change the output, however the asymptotic expansion requires z to lie *
     *  in the right part of the plane.                                       */
    if (z_real < 0.0)
        z = rssringoccs_CDouble_Multiply_Real(-1.0, z);

    /*  The magnitude of the real part determines if we can use the           *
     *  asymptotic expansion or not. Large real parts result in float         *
     *  overflow, giving us infinity.                                         */
    abs_z_real = rssringoccs_Double_Abs(z_real);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (abs_z < 16.0)
    {
        /*  The series is in powers of z^2, so use Horner's method with that. */
        arg = rssringoccs_CDouble_Multiply(z, z);
        bessel_I0 = rssringoccs_CDouble_Poly_Real_Coeffs(
            __taylor_double, __taylor_double_deg, arg
        );

    }

    /*  For larger values, use the asymptotic expansion. MAX_DOUBLE_BASE_E is *
     *  defined in rss_ringoccs_math.h, it's the largest value for a double   *
     *  which won't return INFINITY for exp(x).                               */
    else if (abs_z_real < rssringoccs_Max_Double_Base_E)
    {
        /*  The asymptotic expansion is in terms of 1/z.                      */
        arg = rssringoccs_CDouble_Reciprocal(z);
        bessel_I0 = rssringoccs_CDouble_Poly_Real_Coeffs(
            __asym_double, __asym_double_deg, arg
        );

        /*  Multiply by the coefficient factor and return.                    */
        arg = rssringoccs_CDouble_Multiply_Real(rssringoccs_Two_Pi, z);
        exp_arg = rssringoccs_CDouble_Exp(z);
        sqrt_arg = rssringoccs_CDouble_Sqrt(arg);
        arg = rssringoccs_CDouble_Divide(exp_arg, sqrt_arg);
        bessel_I0 = rssringoccs_CDouble_Multiply(bessel_I0, arg);
    }

    /*  For very large inputs, return INFINITY. INFINITY is standard in C99,  *
     *  and provided in rss_ringoccs_math.h otherwise.                        */
    else
        bessel_I0 = rssringoccs_CDouble_Rect(rssringoccs_Infinity,
                                                   rssringoccs_Infinity);

    return bessel_I0;
}
