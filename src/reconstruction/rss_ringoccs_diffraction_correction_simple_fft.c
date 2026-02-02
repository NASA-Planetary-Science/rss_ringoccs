#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_fft.h>
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stdlib.h>

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionSimpleFFT                                        *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using an FFT about the center of the    *
 *      data. This is the fastest method, but assumes the geometry about the  *
 *      midpoint is an accurate representation of the entire occultation.     *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary     *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This method is the fast, but least accurate. It is very accurate  *
 *          near the midpoint, but assumes the geometry of this point is a    *
 *          fair representative of all of the geometry. The further one gets  *
 *          from the midpoint, the less accurate this is.                     *
 *      2.) This function uses FFTW, which is a NON-STANDARD C Library. Since *
 *          there is no FFT routine in the standard C library, FFTW has       *
 *          somewhat become the de facto standard.                            *
 ******************************************************************************/
void rssringoccs_Diffraction_Correction_SimpleFFT(rssringoccs_TAUObj *tau)
{
    /*  Variables for indexing. nw_pts is the number of points in the window. */
    size_t i, nw_pts, center, data_size, shift, start;

    /*  Some variables needed for reconstruction.                             */
    double w_init, w_thresh;
    tmpl_ComplexDouble *ker;
    tmpl_ComplexDouble *fft_in;
    tmpl_ComplexDouble *fft_out;
    tmpl_ComplexDouble *T_out;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Core_Data(tau);

    /* Check to ensure you have enough data to the left.                      */
    rssringoccs_Tau_Check_Data_Range(tau);

    if (tau->error_occurred)
        return;

    /*  Variable for the center of the data set.                              */
    center = tau->start + tau->n_used / 2;

    /*  Compute some more variables.                                          */
    w_init = tau->w_km_vals[center];
    nw_pts = (TMPL_CAST(w_init / tau->dx_km, size_t) << 1) + 1;
    start = tau->start - (nw_pts - 1) / 2;

    /*  Number of points in the data set, the start/end point and array size. */
    data_size = tau->n_used + nw_pts - 1;

    /* Variables for shifting and keeping track of indexing.                  */
    shift = data_size / 2;

    /*  Allocate memory for the Fresnel kernel and other variables.           */
    ker = malloc(sizeof(*ker) * data_size);
    fft_in = malloc(sizeof(*fft_in) * data_size);

    /*  We can reuse the memory allocated for the in variables for the out    *
     *  variables. This saves us calls to malloc and wasteful memory use.     */
    fft_out = ker;
    T_out = fft_in;

    w_thresh = 0.5*tau->w_km_vals[center];

    /*  Compute the windowing function and Psi.                               */
    for (i = 0; i < data_size; ++i)
    {
        const size_t offset = start + i ;
        const double x = tau->rho_km_vals[offset] - tau->rho_km_vals[center];

        if (tmpl_Double_Abs(x) <= w_thresh)
        {
            const double taper = tau->window_func(x, tau->w_km_vals[center]);

            ker[i] = rssringoccs_Fresnel_Kernel(tau, center, offset);
            tmpl_CDouble_MultiplyBy_Real(&ker[i], taper * tau->dx_km);
        }

        else
            ker[i] = tmpl_CDouble_Zero;
    }

    tmpl_CDouble_FFT(ker, fft_out, data_size);
    tmpl_CDouble_FFT(tau->T_in + start, fft_in, data_size);

    for (i = 0; i < data_size; ++i)
        tmpl_CDouble_MultiplyBy(&fft_out[i], &fft_in[i]);

    tmpl_CDouble_IFFT(fft_out, T_out, data_size);

    for(i = 0; i < tau->n_used; ++i)
    {
        const size_t i_shift = ((nw_pts-1)/2 + i + shift) % (data_size);
        tau->T_out[tau->start + i] = T_out[i_shift];
    }

    /*  Free variables allocated by malloc.                                   */
    free(ker);
    free(fft_in);
}
