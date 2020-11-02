#include <stdio.h>
#include <rss_ringoccs_complex.h>

int main(void)
{
    double real, imag;
    rssringoccs_ComplexDouble z0, z1, erf_z;

    z0 = rssringoccs_Complex_Rect(1.0, 1.0);
    z1 = rssringoccs_Complex_Rect(0.0, 3.1);

    erf_z = rssringoccs_Complex_Erf(z0);
    real = rssringoccs_Complex_Real_Part(erf_z);
    imag = rssringoccs_Complex_Imag_Part(erf_z);
    printf("erf(1+1i) = %f + i%f\n", real, imag);

    erf_z = rssringoccs_Complex_Erf(z1);
    real = rssringoccs_Complex_Real_Part(erf_z);
    imag = rssringoccs_Complex_Imag_Part(erf_z);
    printf("erf(1+1i) = %f + i%f\n", real, imag);

    return 0;
}
