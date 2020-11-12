#include <rss_ringoccs/src/complex/rss_ringoccs_complex.h>
#include <stdio.h>

int main(void)
{
    rssringoccs_ComplexDouble z0, z1, z2;
    rssringoccs_ComplexDouble exp_z0, exp_z1, exp_z2;
    double re, im, exp_re, exp_im;

    z0 = rssringoccs_Complex_One;
    z1 = rssringoccs_Complex_Rect(0.0, 3.1415926);
    z2 = rssringoccs_Complex_Rect(1.0, 1.0);

    exp_z0 = rssringoccs_Complex_Exp(z0);
    exp_z1 = rssringoccs_Complex_Exp(z1);
    exp_z2 = rssringoccs_Complex_Exp(z2);

    re = rssringoccs_Complex_Real_Part(z0);
    im = rssringoccs_Complex_Imag_Part(z0);
    exp_re = rssringoccs_Complex_Real_Part(exp_z0);
    exp_im = rssringoccs_Complex_Imag_Part(exp_z0);
    printf("exp(%f + i%f) = %f + i%f\n", re, im, exp_re, exp_im);

    re = rssringoccs_Complex_Real_Part(z1);
    im = rssringoccs_Complex_Imag_Part(z1);
    exp_re = rssringoccs_Complex_Real_Part(exp_z1);
    exp_im = rssringoccs_Complex_Imag_Part(exp_z1);
    printf("exp(%f + i%f) = %f + i%f\n", re, im, exp_re, exp_im);

    re = rssringoccs_Complex_Real_Part(z2);
    im = rssringoccs_Complex_Imag_Part(z2);
    exp_re = rssringoccs_Complex_Real_Part(exp_z2);
    exp_im = rssringoccs_Complex_Imag_Part(exp_z2);
    printf("exp(%f + i%f) = %f + i%f\n", re, im, exp_re, exp_im);

    return 0;
}
