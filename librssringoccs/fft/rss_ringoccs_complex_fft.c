#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fft.h>

rssringoccs_ComplexDouble *
rssringoccs_Complex_FFT(rssringoccs_ComplexDouble *in, unsigned long N,
                        rssringoccs_Bool inverse)
{
    rssringoccs_ComplexDouble *out;
    if ((N > 0) && ((N & (N-1)) == 0))
        out = rssringoccs_Complex_FFT_Cooley_Tukey(in, N, inverse);
    else
        out = rssringoccs_Complex_FFT_Bluestein_Chirp_Z(in, N, inverse);

    return out;
}
