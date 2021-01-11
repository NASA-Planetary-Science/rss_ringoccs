#include <stdlib.h>
#include <stdio.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fft.h>

rssringoccs_ComplexDouble *
rssringoccs_Complex_FFT_Bluestein_Chirp_Z(rssringoccs_ComplexDouble *in,
                                          unsigned long N,
                                          rssringoccs_Bool inverse)
{
    /*  The chirp factors range from -(N-1) to N-1, inclusive.                */
    unsigned long chirp_size;

    /*  We're going to use the radix-2 FFT to compute the FFT of an arbitrary *
     *  array. We'll need the next highest power of 2 greater than N. Use     *
     *  bitwise operations to do this.                                        */
    unsigned long N_pow_2;

    /*  Variables for indexing.                                               */
    unsigned long n, m;

    /*  We'll need some variables declared. These will be our arrays.         */
    rssringoccs_ComplexDouble *chirp;
    rssringoccs_ComplexDouble *rcpr_chirp;
    rssringoccs_ComplexDouble *x_in;
    rssringoccs_ComplexDouble *fft_x_in;
    rssringoccs_ComplexDouble *fft_rcpr_chirp;
    rssringoccs_ComplexDouble *out;

    /*  And scale factors for the exponential.                                */
    double chirp_factor;
    double m2_chirp_factor;
    double inverse_factor;

    chirp_size = N+N-1;

    /*  Now, to get the highest power of two greater than N, think of how you *
     *  would do it for the highest power of 10. You would simply count off   *
     *  the number of digits. For example, if N = 1436, then there are four   *
     *  digits and the next highest power of 10 larger than N is 10^4 = 10000.*
     *  Do the same thing, but in binary! We do this as follows:              */

    /*  First, set N_pow_2 to zero so it has been given a starting value.     */
    N_pow_2 = 0;
    n = 0;

    /*  Now we count the number of digits in the binary representation of N.  */
    while(N_pow_2 < chirp_size)
    {
        /*  Incremenent i.                                                    */
        ++n;

        /*  We're going to set N_pow_2 to 1000...000 where there are i zeros. */
        N_pow_2 = 1ul << n;
    }

    /*  Allocate memory for x_in and chirp, which will be a power of two in   *
     *  size. Per C90 guidelines, we do not cast malloc since void pointers   *
     *  safely promoted without the need for type casting.                    */
    chirp      = malloc(sizeof(*chirp)      * chirp_size);
    rcpr_chirp = malloc(sizeof(*rcpr_chirp) * N_pow_2);
    x_in       = malloc(sizeof(*x_in)       * N_pow_2);

    if (inverse)
        chirp_factor = rssringoccs_One_Pi/(double)N;
    else
        chirp_factor = -rssringoccs_One_Pi/(double)N;

    /*  Set the values for the "chirp" factor, which is simply the complex    *
     *  exponential of (k^2 / 2) * (+/- 2 pi i / N). The +/- depends on       *
     *  whether or not an inverse computation is being performed.             */
    for (n=0; n<chirp_size; ++n)
    {
        m = n+1-N;
        m2_chirp_factor = (double)(m*m)*chirp_factor;
        chirp[n] = rssringoccs_CDouble_Polar(1.0, m2_chirp_factor);
        rcpr_chirp[n] = rssringoccs_CDouble_Reciprocal(chirp[n]);
    }

    /*  Now pad the rest of chirp with zeros so that it is a power of two.    */
    for (n=chirp_size; n<N_pow_2; ++n)
        rcpr_chirp[n] = rssringoccs_CDouble_Zero;

    /*  Set the x_in array to in times chirp, and then pad with zero.         */
    for (n=0; n<N; ++n)
        x_in[n] = rssringoccs_CDouble_Multiply(chirp[n+N-1], in[n]);

    /*  Now pad the rest with zeros.                                          */
    for (n=N; n<N_pow_2; ++n)
        x_in[n] = rssringoccs_CDouble_Zero;

    /*  Lastly, we need to compute the forward FFTs of x_in and chirp, and    *
     *  then compute the inverse fourier transform of the product. We'll need *
     *  to allocate memory for these two.                                     */
    fft_x_in = rssringoccs_Complex_FFT_Cooley_Tukey(x_in, N_pow_2,
                                                    rssringoccs_False);

    if (fft_x_in == NULL)
    {
        puts("fft_x_in 1 failed.");
        free(x_in);
        free(chirp);
        free(rcpr_chirp);
        return NULL;
    }

    fft_rcpr_chirp = rssringoccs_Complex_FFT_Cooley_Tukey(
        rcpr_chirp, N_pow_2, rssringoccs_False
    );

    if (fft_rcpr_chirp == NULL)
    {
        puts("fft_rcpr_chirp 1 failed.");
        free(x_in);
        free(chirp);
        free(rcpr_chirp);
        free(fft_x_in);
        return NULL;
    }

    for (n=0; n<N_pow_2; ++n)
        x_in[n] = rssringoccs_CDouble_Multiply(fft_x_in[n], fft_rcpr_chirp[n]);

    free(fft_x_in);
    fft_x_in = rssringoccs_Complex_FFT_Cooley_Tukey(x_in, N_pow_2,
                                                    rssringoccs_True);

    if (fft_x_in == NULL)
    {
        puts("fft_x_in 2 failed.");
        free(x_in);
        free(chirp);
        free(fft_rcpr_chirp);
        free(fft_x_in);
        return NULL;
    }

    out = malloc(sizeof(*out) * N);

    for(n=0; n<N; ++n)
    {
        m = n+N-1;
        out[n] = rssringoccs_CDouble_Multiply(fft_x_in[m], chirp[m]);
    }

    if (inverse)
    {
        inverse_factor = 1.0/(double)N;
        for (n=0; n<N; ++n)
            out[n] = rssringoccs_CDouble_Multiply_Real(inverse_factor, out[n]);
    }

    /*  Don't forget to free everything!!!                                    */
    free(x_in);
    free(chirp);
    free(rcpr_chirp);
    free(fft_x_in);
    free(fft_rcpr_chirp);

    return out;
}
