#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>

rssringoccs_ComplexDouble *
rssringoccs_FFT_Cooley_Tukey_ComplexDouble(rssringoccs_ComplexDouble *in,
                                           long N, rssringoccs_Bool inverse)
{
    /*  We'll need several variables for indexing the pointers.               */
    long k, m, n, skip;

    /*  Boolean for determining if we're on an even iteration.                */
    rssringoccs_Bool evenIteration = N & 0x55555555;

    /*  Declare several pointers for performing the Cooley-Tukey algorithm.   */
    rssringoccs_ComplexDouble *out;
    rssringoccs_ComplexDouble *E, *D;
    rssringoccs_ComplexDouble *twiddles;
    rssringoccs_ComplexDouble *scratch;
    rssringoccs_ComplexDouble *Xp;
    rssringoccs_ComplexDouble *Xp2;
    rssringoccs_ComplexDouble *Xstart;

    /*  And some variables for actually computing the FFT.                    */
    rssringoccs_ComplexDouble t, d;
    double factor;

    /*  This method assume N is a power of two. If not, return failure.       */
    if (!((N > 0) && ((N & (N-1)) == 0)))
        return NULL;

    out = malloc(sizeof(*out) * N);

    /*  If N is 1, simply return the output. That is, the FFT of a point is   *
     *  just that point.                                                      */
    if (N == 1)
    {
    	out[0] = in[0];
    	return out;
    }

    /*  The "twiddle" factors are just the complex exponentials that occur    *
     *  inside the discrete Fourier transform. Allocate memory for this and   *
     *  the "scratch" factor. Per C99 recommendations we do not cast malloc.  */ 
    twiddles = malloc(sizeof(*twiddles) * N);
    scratch  = malloc(sizeof(*scratch)  * N);

    /*  If we are performing an inverse Fourier transform, the factor inside  *
     *  the exponential is 2 pi / N. Forward transform is minus this.         */
    if (inverse)
        factor = 2.0 * M_PI/N;
    else
        factor = -2.0 * M_PI/N;

    /*  Compute the "twiddle" factors. No idea why it's called this.          */
    for (k = 0; k<N; ++k)
        twiddles[k] = rssringoccs_Complex_Polar(1.0, k * factor);

    /*  Set "E" pointer to point to the initial address of the "in" pointer.  */
    E = in;

    /*  The actual Cooley-Tukey algorithm is recursive. We can save a lot of  *
     *  overhall, and memory usage, by unraveling this a bit into nested      *
     *  for-loops. This is also significantly faster than pure recursion.     */
    for (n = 1; n < N; n *= 2)
    {
        if (evenIteration)
            Xstart = scratch;
        else
            Xstart = out;

        skip = N/(2 * n);

        /* Each of D and E is of length n, and each element of each D and E   *
         * is separated by 2*skip. The Es begin at E[0] to E[skip - 1] and    *
         * the Ds begin at E[skip] to E[2*skip - 1]                           */
        Xp = Xstart;
        Xp2 = Xstart + N/2;
        for (k=0; k<n; ++k)
        {
        	t = twiddles[k * skip];
            for (m=0; m<skip; ++m)
            {
                /*  Set the D pointer to the desired address.                 */
            	D = E + skip;

            	/* twiddle *D to get dre and dim                              */
            	d = rssringoccs_Complex_Multiply(t, *D);
                *Xp  = rssringoccs_Complex_Add(*E, d);
                *Xp2 = rssringoccs_Complex_Subtract(*E, d);
                ++Xp;
                ++Xp2;
                ++E;
            }
            E += skip;
        }
        E = Xstart;

        /*  The next iteration is the opposite of what evenIteration is now.  */
        evenIteration = !evenIteration;
    }

    /*  Free your allocated memory!                                           */
    free(twiddles);
    free(scratch);

    /*  The inverse Fourier transform has a 1/N factor in front of the sum.   */
    if (inverse)
    {
        factor = 1.0/N;
        for (k=0; k<N; ++k)
            out[k] = rssringoccs_Complex_Scale(factor, out[k]);
    }

    return out;
}

rssringoccs_ComplexDouble *
rssringoccs_FFT_Bluestein_Chirp_Z_ComplexDouble(rssringoccs_ComplexDouble *in,
                                                long N,
                                                rssringoccs_Bool inverse)
{
    /*  The chirp factors range from -(N-1) to N-1, inclusive.                */
    long chirp_size;

    /*  We're going to use the radix-2 FFT to compute the FFT of an arbitrary *
     *  array. We'll need the next highest power of 2 greater than N. Use     *
     *  bitwise operations to do this.                                        */
    long N_pow_2;

    /*  Variables for indexing.                                               */
    long n, m;

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

    /*  N should be a positive integer.                                       */
    if (N<0)
        return NULL;

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
        N_pow_2 = 1 << n;
    }

    /*  Allocate memory for x_in and chirp, which will be a power of two in   *
     *  size. Per C90 guidelines, we do not cast malloc since void pointers   *
     *  safely promoted without the need for type casting.                    */
    chirp      = malloc(sizeof(*chirp)      * chirp_size);
    rcpr_chirp = malloc(sizeof(*rcpr_chirp) * N_pow_2);
    x_in       = malloc(sizeof(*x_in)       * N_pow_2);

    if (inverse)
        chirp_factor = M_PI/N;
    else
        chirp_factor = -M_PI/N;

    /*  Set the values for the "chirp" factor, which is simply the complex    *
     *  exponential of (k^2 / 2) * (+/- 2 pi i / N). The +/- depends on       *
     *  whether or not an inverse computation is being performed.             */
    for (n=0; n<chirp_size; ++n)
    {
        m = n+1-N;
        m2_chirp_factor = m*m*chirp_factor;
        chirp[n] = rssringoccs_Complex_Polar(1.0, m2_chirp_factor);
        rcpr_chirp[n] = rssringoccs_Complex_Reciprocal(chirp[n]);
    }

    /*  Now pad the rest of chirp with zeros so that it is a power of two.    */
    for (n=chirp_size; n<N_pow_2; ++n)
        rcpr_chirp[n] = rssringoccs_Complex_Zero;

    /*  Set the x_in array to in times chirp, and then pad with zero.         */
    for (n=0; n<N; ++n)
        x_in[n] = rssringoccs_Complex_Multiply(chirp[n+N-1], in[n]);

    /*  Now pad the rest with zeros.                                          */
    for (n=N; n<N_pow_2; ++n)
        x_in[n] = rssringoccs_Complex_Zero;

    /*  Lastly, we need to compute the forward FFTs of x_in and chirp, and    *
     *  then compute the inverse fourier transform of the product. We'll need *
     *  to allocate memory for these two.                                     */
    fft_x_in = rssringoccs_FFT_Cooley_Tukey_ComplexDouble(x_in, N_pow_2,
                                                          rssringoccs_False);

    if (!fft_x_in)
    {
        free(x_in);
        free(chirp);
        return NULL;
    }

    fft_rcpr_chirp = rssringoccs_FFT_Cooley_Tukey_ComplexDouble(
        rcpr_chirp, N_pow_2, rssringoccs_False
    );

    if (!fft_x_in)
    {
        free(x_in);
        free(chirp);
        return NULL;
    }
    
    for (n=0; n<N_pow_2; ++n)
        x_in[n] = rssringoccs_Complex_Multiply(fft_x_in[n], fft_rcpr_chirp[n]);

    fft_x_in = rssringoccs_FFT_Cooley_Tukey_ComplexDouble(x_in, N_pow_2,
                                                          rssringoccs_True);

    if (!fft_x_in)
    {
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
        out[n] = rssringoccs_Complex_Multiply(fft_x_in[m], chirp[m]);
    }

    if (inverse)
    {
        inverse_factor = 1.0/N;
        for (n=0; n<N; ++n)
            out[n] = rssringoccs_Complex_Scale(inverse_factor, out[n]);
    }

    /*  Don't forget to free everything!!!                                    */
    free(x_in);
    free(chirp);
    free(rcpr_chirp);
    free(fft_x_in);
    free(fft_rcpr_chirp);

    return out;
}
