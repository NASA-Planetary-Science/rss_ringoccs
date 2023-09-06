#include <stdlib.h>
#include <math.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_fft.h>
#include <libtmpl/include/tmpl_cyl_fresnel_optics.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionSimpleFFT                                        *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using an FFT about the center of the    *
 *      data. This is the fastest method, but assumes the geometry about the  *
 *      midpoint is an accurate representation of the entire occultation.
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary    *
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
    size_t i, nw_pts, center, data_size, shift, current_point, i_shift;

    /*  Some variables needed for reconstruction.                             */
    double w_init, psi, phi, window_func_x, factor, rcpr_F;
    double w_thresh, arg_norm, D;
    tmpl_ComplexDouble *ker;
    tmpl_ComplexDouble *fft_ker;
    tmpl_ComplexDouble *fft_in;
    tmpl_ComplexDouble *fft_out;
    tmpl_ComplexDouble *T_in;
    tmpl_ComplexDouble *T_out;
    tmpl_ComplexDouble arg;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Data(tau);
    if (tau->error_occurred)
        return;

    /* Check to ensure you have enough data to the left.                      */
    rssringoccs_Tau_Check_Data_Range(tau);
    if (tau->error_occurred)
        return;

    /*  Variable for the center of the data set.                              */
    center = tau->start + tau->n_used/2;

    /*  Compute some more variables.                                          */
    w_init = tau->w_km_vals[center];
    nw_pts = (size_t)(w_init / (tau->dx_km * 2.0));

    /*  Number of points in the data set, the start/end point and array size. */
    data_size = tau->n_used + 2*nw_pts + 1;

    /* Variables for shifting and keeping track of indexing.                  */
    shift = data_size / 2;

    /*  Scale factor for the FFT.                                             */
    factor = 0.5*tau->dx_km;

    /*  Allocate memory for the Fresnel kernel and other variables.           */
    ker     = malloc(sizeof(*ker)     * data_size);
    T_in    = malloc(sizeof(*T_in)    * data_size);

    /*  We can reuse the memory allocated for the in variables for the out    *
     *  variables. This saves four calls to malloc and redundant memory use.  */
    fft_in = T_in;
    fft_ker = ker;
    fft_out = fft_in;
    T_out = T_in;

    w_thresh = 0.5*tau->w_km_vals[center];

    /*  Compute the windowing function and Psi.                               */
    for (i=0; i < data_size; ++i)
    {
        current_point = tau->start + i - nw_pts;
        window_func_x = tau->rho_km_vals[center] -
                        tau->rho_km_vals[current_point];

        if (fabs(window_func_x) <= w_thresh)
        {
            phi = tmpl_Double_Stationary_Cyl_Fresnel_Psi_D_Newton(
                tau->k_vals[center],
                tau->rho_km_vals[center],
                tau->rho_km_vals[current_point],
                tau->phi_deg_vals[current_point],
                tau->phi_deg_vals[current_point],
                tau->B_deg_vals[center],
                tau->rx_km_vals[center],
                tau->ry_km_vals[center],
                tau->rz_km_vals[center],
                tau->EPS,
                tau->toler
            );

            D = tmpl_Double_Cyl_Fresnel_Observer_Distance(
                tau->rho_km_vals[current_point],   /* Ring radius. */
                phi,                               /* Stationary azimuth. */
                tau->rx_km_vals[center],           /* Cassini x coordinate. */
                tau->ry_km_vals[center],           /* Cassini y coordinate. */
                tau->rz_km_vals[center]            /* Cassini z coordinate. */
            );

            psi = -tmpl_Double_Cyl_Fresnel_Psi(
                tau->k_vals[center],
                tau->rho_km_vals[center],
                tau->rho_km_vals[current_point],
                phi,
                tau->phi_deg_vals[current_point],
                tau->B_deg_vals[center],
                D
            );

            /*  If forward tranform is set, negate the Fresnel kernel.        */
            if (tau->use_fwd)
                psi *= -1.0;

            arg_norm = tau->window_func(window_func_x, tau->w_km_vals[center]);
            ker[i] = tmpl_CDouble_Polar(arg_norm, psi);
        }
        else
            ker[i] = tmpl_CDouble_Zero;


        T_in[i] = tau->T_in[current_point];
    }

    tmpl_CDouble_FFT(ker, fft_ker, data_size);
    tmpl_CDouble_FFT(T_in, fft_in, data_size);

    for (i = 0; i < data_size; ++i)
        fft_out[i] = tmpl_CDouble_Multiply(fft_ker[i], fft_in[i]);

    tmpl_CDouble_IFFT(fft_out, T_out, data_size);

    for(i = 0; i < tau->n_used; ++i)
    {
        i_shift = (nw_pts + i + shift) % (data_size);
        rcpr_F = 1.0/tau->F_km_vals[tau->start + i];
        arg = tmpl_CDouble_Rect(factor*rcpr_F, factor*rcpr_F);
        tau->T_out[tau->start + i] = tmpl_CDouble_Multiply(arg, T_out[i_shift]);
    }

    /*  Free variables allocated by malloc.                                   */
    free(ker);
    free(T_in);
}
