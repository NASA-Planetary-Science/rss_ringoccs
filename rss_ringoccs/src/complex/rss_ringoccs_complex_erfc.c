/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

static double __taylor[51] =
{
    -1.1283791670955126,
    0.3761263890318375,
    -0.11283791670955126,
    0.02686617064513125,
    -0.005223977625442187,
    0.0008548327023450852,
    -0.00012055332981789664,
    0.000014925650358406252,
    -1.6462114365889246e-6,
    1.6365844691234924e-7,
    -1.4807192815879218e-8,
    1.2290555301717926e-9,
    -9.422759064650411e-11,
    6.7113668551641105e-12,
    -4.463224263286477e-13,
    2.7835162072109215e-14,
    -1.6342614095367152e-15,
    9.063970842808672e-17,
    -4.763348040515068e-18,
    2.3784598852774293e-19,
    -1.131218725924631e-20,
    5.136209054585811e-22,
    -2.230878680274645e-23,
    9.28672901131906e-25,
    -3.71153285316323e-26,
    1.4263930180784176e-27,
    -5.279103332510835e-29,
    1.8841244217042036e-30,
    -6.492909974544561e-32,
    2.1630383901171243e-33,
    -6.973730328792914e-35,
    2.1781748594796097e-36,
    -6.597356545539202e-38,
    1.9395213725013484e-39,
    -5.539127534424142e-41,
    1.538027363683162e-42,
    -4.1552489658106737e-44,
    1.0930925207357809e-45,
    -2.8018434400267797e-47,
    7.002335114640116e-49,
    -1.7073594878289173e-50,
    4.0639470618319806e-52,
    -9.448392328628974e-54,
    2.1467878854142287e-55,
    -4.769421502324767e-57,
    1.0365775670498273e-58,
    -2.2049686442621383e-60,
    4.59265585479012e-62,
    -9.3707539992496e-64,
    1.8737644566629794e-65,
    -3.673320419992772e-67
};

static double __asym[4] =
{
    1.12838, -1.12838, 3.38514, -16.9257
};

rssringoccs_ComplexDouble
rssringoccs_Complex_Erfc(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires this at the top.            */
    rssringoccs_ComplexDouble erfcz, arg, exp_arg;
    double abs_z, abs_z_real, z_real;
    int deg;

    /*  Extract the real part of z.                                           */
    z_real = rssringoccs_Complex_Real_Part(z);

    /*  Compute the magnitude of z and the absolute value of its real part.   */
    abs_z = rssringoccs_Complex_Abs(z);

    /*  If z is negative, flip z. That is, negative z so that it lies in the  *
     *  right-half of the complex plane. I0 is even so negation does not      *
     *  change the output, however the asymptotic expansion requires z to lie *
     *  in the right part of the plane.                                       */
    if (z_real < 0.0)
        z = rssringoccs_Complex_Scale(-1.0, z);

    /*  The magnitude of the real part determines if we can use the           *
     *  asymptotic expansion or not. Large real parts result in float         *
     *  overflow, giving us infinity.                                         */
    abs_z_real = rssringoccs_Abs_Double(z_real);

    /*  For small arguments, use a Taylor series to approximate I_0.          */
    if (abs_z < 5.0)
    {
        arg = rssringoccs_Complex_Multiply(z, z);
        if (abs_z < 1.0)
            deg = 10;
        else if (abs_z < 2.0)
            deg = 20;
        else if (abs_z < 3.0)
            deg = 30;
        else if (abs_z < 4.0)
            deg = 40;
        else
            deg = 50;

        erfcz = rssringoccs_Complex_Poly_Real_Coeffs(__taylor, deg, arg);
        erfcz = rssringoccs_Complex_Multiply(z, erfcz);
    }

    /*  For larger values, use the asymptotic expansion. MAX_DOUBLE_BASE_E is *
     *  defined in rss_ringoccs_math.h, it's the largest value for a double   *
     *  which won't return INFINITY for exp(x).                               */
    else
    {
        /*  The asymptotic expansion is in terms of 1/z.                      */
        arg = rssringoccs_Complex_Reciprocal(z);
        erfcz = rssringoccs_Complex_Poly_Real_Coeffs(__asym, 3, arg);
        erfcz = rssringoccs_Complex_Multiply(arg, erfcz);

        /*  Multiply by the coefficient factor and return.                    */
        arg = rssringoccs_Complex_Multiply(z, z);
        arg = rssringoccs_Complex_Scale(-1.0, arg);
        exp_arg = rssringoccs_Complex_Exp(arg);
        erfcz = rssringoccs_Complex_Multiply(exp_arg, erfcz);
    }

    return erfcz;
}
