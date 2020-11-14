/*  Header file which contains aliases for the function in the standard C     *
 *  library math.h. This allows compatibility of C89 and C99 math.h headers.  */
#include <rss_ringoccs/include/rss_ringoccs_math.h>

/*  Where the prototypes are declared and where complex types are defined.    */
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

/*  If C99 complex.h was included before rss_ringoccs_complex.h then we'll    *
 *  use the built-in +, *, -, and / for complex arithmetic.                   */
#if _RSS_RINGOCCS_USING_COMPLEX_H_ == 1

const rssringoccs_ComplexDouble rssringoccs_Imaginary_Unit = _Complex_I;
const rssringoccs_ComplexDouble rssringoccs_Complex_Zero = 0.0;
const rssringoccs_ComplexDouble rssringoccs_Complex_One = 1.0;

/*  Function for creating a new complex number from two real numbers.         */
rssringoccs_ComplexDouble rssringoccs_Complex_Rect(double x, double y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z;

    /*  We simply need the complex number x+iy. This can be computed with the *
     *  built-in + and * in C99.                                              */
    z = x + rssringoccs_Imaginary_Unit*y;
    return z;
}

rssringoccs_ComplexDouble rssringoccs_Complex_Polar(double r, double theta)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real, imag;
    rssringoccs_ComplexDouble z;

    /*  We use Euler's formula to compute z = r*exp(i*theta) and write this   *
     *  as z = r*cos(theta) + i*r*sin(theta). We compute the real and         *
     *  imaginary parts and then combine them.                                */
    real = r*rssringoccs_Cos_Double(theta);
    imag = r*rssringoccs_Sin_Double(theta);
    z = real + rssringoccs_Imaginary_Unit*imag;
    return z;
}

/*  Function for adding two complex numbers.                                  */
rssringoccs_ComplexDouble rssringoccs_Complex_Add(rssringoccs_ComplexDouble z0,
                                                  rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sum;

    /*  + works on complex numbers in C99, so simply use this and return.     */
    sum = z0+z1;
    return sum;
}

/*  Function for subtracting two complex numbers.                             */
rssringoccs_ComplexDouble
rssringoccs_Complex_Subtract(rssringoccs_ComplexDouble z0,
                             rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble diff;

    /*  - works on complex numbers in C99, so simply use this and return.     */
    diff = z0-z1;
    return diff;
}

/*  Function for multiplying two complex numbers.                             */
rssringoccs_ComplexDouble
rssringoccs_Complex_Multiply(rssringoccs_ComplexDouble z0,
                             rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble prod;

    /*  * works on complex numbers in C99, so simply use this and return.     */
    prod = z0*z1;
    return prod;
}

/*  Function for multiplying a complex number by a real numbers.              */
rssringoccs_ComplexDouble
rssringoccs_Complex_Scale(double x, rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble scale;

    /*  C99 allows one to mix * with real and complex data types, so do this. */
    scale = x*z;
    return scale;
}

/*  Function for computing 1/z for non-zero z.                                */
rssringoccs_ComplexDouble
rssringoccs_Complex_Reciprocal(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble rcpr_z;

    /*  C99 allows use of the / operator with complex numbers, so use this.   */
    rcpr_z = 1.0/z;
    return rcpr_z;
}

/*  Function for dividing a complex number by another.                        */
rssringoccs_ComplexDouble
rssringoccs_Complex_Divide(rssringoccs_ComplexDouble z0,
                           rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble quotient;

    /*  / works on complex numbers in C99, so simply use this and return.     */
    quotient = z0/z1;
    return quotient;
}

/*  If you do not have complex.h supported, we'll need to define complex      *
 *  arithmetic using the rssringoccs_ComplexDouble data type.                 */
#else

/*  Useful constants.                                                         */
const rssringoccs_ComplexDouble rssringoccs_Imaginary_Unit = {{0.0, 1.0}};
const rssringoccs_ComplexDouble rssringoccs_Complex_Zero   = {{0.0, 0.0}};
const rssringoccs_ComplexDouble rssringoccs_Complex_One    = {{1.0, 0.0}};

/*  In C99 you can simply do double _Complex z = x + _Complex_I*y since       *
 *  complex variables are primitive data types, but in C89 we need to create  *
 *  a struct for them. Structs can't be added, so we need a function for      *
 *  creating a complex number from two doubles.                               */

/*  Function for creating a complex number from its real and imaginary parts. */
rssringoccs_ComplexDouble rssringoccs_Complex_Rect(double x, double y)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z;

    /*  Simply set the array dat inside z to {x, y} and return.               */
    z.dat[0] = x;
    z.dat[1] = y;
    return z;
}

/*  Returns a complex number given a polar representation (r, theta).         */
rssringoccs_ComplexDouble rssringoccs_Complex_Polar(double r, double theta)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble z;
    double real, imag;

    /*  Use Euler's formula for the polar representation of a complex number. */
    real = r * rssringoccs_Cos_Double(theta);
    imag = r * rssringoccs_Sin_Double(theta);

    /*  Use rssringoccs_ComplexRect to compute the complex number and return. */
    z = rssringoccs_Complex_Rect(real, imag);
    return z;
}

/*  This function is equivalent to the creal function in complex.h (C99).     */
double rssringoccs_Complex_Real_Part(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real;

    /*  The real component is stored as the first entry in the dat array      *
     *  contained in a rssringoccs_ComplexDouble struct. Return this.         */
    real = z.dat[0];
    return real;
}

/*  This function is equivalent to the cimag function in complex.h (C99).     */
double rssringoccs_Complex_Imag_Part(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double imag;

    /*  The imaginary component is stored as the second entry in the dat      *
     *  array contained in a rssringoccs_ComplexDouble struct. Return this.   */
    imag = z.dat[1];
    return imag;
}

rssringoccs_ComplexDouble
rssringoccs_Complex_Conjugate(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real, imag;
    rssringoccs_ComplexDouble conjz;

    /*  Extract the values from the complex number.                           */
    real = rssringoccs_Complex_Real_Part(z);
    imag = rssringoccs_Complex_Imag_Part(z);

    /*  The complex conjugate of x+iy is just x-iy, compute this.             */
    conjz = rssringoccs_Complex_Rect(real, -imag);
    return conjz;
}

/*  This function is equivalent to the cabs function in complex.h (C99). It   *
 *  computes the absolute value of a complex number z.                        */
double rssringoccs_Complex_Abs(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real, imag, abs_value;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_Complex_Real_Part(z);
    imag = rssringoccs_Complex_Imag_Part(z);

    /*  The absolute value is just sqrt(x^2 + y^2) so compute this.           */
    abs_value = rssringoccs_Sqrt_Double(real*real + imag*imag);
    return abs_value;
}

/*  Function for computing the angle of a complex number. Returns a value in  *
 *  the range (-pi, pi].                                                      */
double rssringoccs_Complex_Argument(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    double real, imag, theta;

    /*  Extract the real and imaginary parts from the input complex number.   */
    real = rssringoccs_Complex_Real_Part(z);
    imag = rssringoccs_Complex_Imag_Part(z);

    /*  Compute the argument using arctan and return.                         */
    theta = rssringoccs_Arctan_Double(imag, real);
    return theta;
}

/*  In C99, since _Complex is a built-in data type, given double _Complex z1  *
 *  and double _Complex z2, you can just do z1 + z2. Structs can't be added,  *
 *  so we need a function for computing the sum of two complex values.        */
rssringoccs_ComplexDouble
rssringoccs_Complex_Add(rssringoccs_ComplexDouble z0,
                        rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble sum;
    double real0, real1;
    double imag0, imag1;
    double sumre, sumim;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_Complex_Real_Part(z0);
    real1 = rssringoccs_Complex_Real_Part(z1);

    imag0 = rssringoccs_Complex_Imag_Part(z0);
    imag1 = rssringoccs_Complex_Imag_Part(z1);

    /*  The sum of two complex numbers simply adds their components.          */
    sumre = real0 + real1;
    sumim = imag0 + imag1;

    /*  Create the output from sumre and sumim and return.                    */
    sum = rssringoccs_Complex_Rect(sumre, sumim);
    return sum;
}

rssringoccs_ComplexDouble
rssringoccs_Complex_Subtract(rssringoccs_ComplexDouble z0,
                             rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble diff;
    double real0, real1;
    double imag0, imag1;
    double diffre, diffim;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_Complex_Real_Part(z0);
    real1 = rssringoccs_Complex_Real_Part(z1);

    imag0 = rssringoccs_Complex_Imag_Part(z0);
    imag1 = rssringoccs_Complex_Imag_Part(z1);

    /*  The difference of two complex numbers simply subtracts components.    */
    diffre = real0 - real1;
    diffim = imag0 - imag1;

    /*  Create the output from diffre and diffim and return.                  */
    diff = rssringoccs_Complex_Rect(diffre, diffim);
    return diff;
}

rssringoccs_ComplexDouble
rssringoccs_Complex_Multiply(rssringoccs_ComplexDouble z0,
                             rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble prod;
    double real0, real1;
    double imag0, imag1;
    double prodre, prodim;

    /*  Extract the real and imaginary parts from the inputs.                 */
    real0 = rssringoccs_Complex_Real_Part(z0);
    real1 = rssringoccs_Complex_Real_Part(z1);

    imag0 = rssringoccs_Complex_Imag_Part(z0);
    imag1 = rssringoccs_Complex_Imag_Part(z1);

    /*  The product uses the distributive law in combination with the fact    *
     *  that i^2 = -1. This gives us the following formulas:                  */
    prodre = real0*real1 - imag0*imag1;
    prodim = real0*imag1 + imag0*real1;

    /*  Create the output from prodre and prodim and return.                  */
    prod = rssringoccs_Complex_Rect(prodre, prodim);
    return prod;
}

rssringoccs_ComplexDouble
rssringoccs_Complex_Scale(double x, rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble scale;
    double real, imag;
    double scalereal, scaleimag;

    /*  Extract the real and imaginary parts from z.                          */
    real = rssringoccs_Complex_Real_Part(z);
    imag = rssringoccs_Complex_Imag_Part(z);

    /*  Scale the components by x and return.                                 */
    scalereal = x*real;
    scaleimag = x*imag;
    scale = rssringoccs_Complex_Rect(scalereal, scaleimag);
    return scale;
}

/*  Function for computing the reciprocal (or inverse) of a complex number.   */
rssringoccs_ComplexDouble
rssringoccs_Complex_Reciprocal(rssringoccs_ComplexDouble z)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble conj_z, rcpr_z;
    double abs_z, rcp_abs_z_sq;

    /*  Compute the conjugate of z and its absolute value.                    */
    conj_z = rssringoccs_Complex_Conjugate(z);
    abs_z = rssringoccs_Complex_Abs(z);

    /*  The inverse of z is conj_z / abs_z^2, so return this.                 */
    rcp_abs_z_sq = 1.0/(abs_z*abs_z);
    rcpr_z = rssringoccs_Complex_Scale(rcp_abs_z_sq, conj_z);
    return rcpr_z;
}

rssringoccs_ComplexDouble
rssringoccs_Complex_Divide(rssringoccs_ComplexDouble z0,
                           rssringoccs_ComplexDouble z1)
{
    /*  Declare necessary variables. C89 requires declarations at the top.    */
    rssringoccs_ComplexDouble div, rcpr_z1;

    /*  We'll compute z0/z1 as z0 * (1/z1), which is the reciprocal of z1.    */
    rcpr_z1 = rssringoccs_Complex_Reciprocal(z1);
    div = rssringoccs_Complex_Multiply(z0, rcpr_z1);
    return div;
}

#endif
/*  End of #if _RSS_RINGOCCS_USING_COMPLEX_H_ == 1                            */
