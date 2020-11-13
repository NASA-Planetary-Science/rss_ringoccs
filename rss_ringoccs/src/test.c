#include <rss_ringoccs/src/complex/rss_ringoccs_complex.h>
#include <stdio.h>

int main(void)
{
    rssringoccs_ComplexDouble z, pow;
    double re, im;

    pow = rssringoccs_Complex_Pow(rssringoccs_Imaginary_Unit,
                                  rssringoccs_Imaginary_Unit);
    re = rssringoccs_Complex_Real_Part(pow);
    im = rssringoccs_Complex_Imag_Part(pow);
    printf("i^i = %f + i%f\n", re, im);

    z = rssringoccs_Complex_Rect(0.7071067811865476,
                                 0.7071067811865476);
    pow = rssringoccs_Complex_Real_Pow(z, 2);
    re = rssringoccs_Complex_Real_Part(pow);
    im = rssringoccs_Complex_Imag_Part(pow);
    printf("((1+i)/sqrt(2))^2 = %f + i%f\n", re, im);
    return 0;
}

