#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_fft.h>

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
    unsigned long i, nw_pts, center, data_size;
    unsigned long shift, current_point, i_shift;

    /*  Some variables needed for reconstruction.                             */
    double w_init, psi, phi, window_func_x, factor, rcpr_F;
    double w_thresh, arg_norm, D, x, y, z, dx, dy;
    rssringoccs_ComplexDouble *ker;
    rssringoccs_ComplexDouble *fft_ker;
    rssringoccs_ComplexDouble *fft_in;
    rssringoccs_ComplexDouble *fft_out;
    rssringoccs_ComplexDouble *T_in;
    rssringoccs_ComplexDouble *T_out;
    rssringoccs_ComplexDouble arg;

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
    nw_pts = (unsigned long)(w_init / (tau->dx_km * 2.0));

    /*  Number of points in the data set, the start/end point and array size. */
    data_size = tau->n_used + 2*nw_pts + 1;

    /* Variables for shifting and keeping track of indexing.                  */
    shift = data_size / 2;

    /*  Scale factor for the FFT.                                             */
    factor = 0.5*tau->dx_km;

    /*  Allocate memory for the Fresnel kernel and other variables.           */
    ker     = malloc(sizeof(*ker)     * data_size);
    fft_out = malloc(sizeof(*fft_out) * data_size);
    T_in    = malloc(sizeof(*T_in)    * data_size);

    w_thresh = 0.5*tau->w_km_vals[center];

    /*  Compute the windowing function and Psi.                               */
    for (i=0; i < data_size; ++i)
    {
        current_point = tau->start + i - nw_pts;
        window_func_x = tau->rho_km_vals[center] -
                        tau->rho_km_vals[current_point];
        if (fabs(window_func_x) <= w_thresh)
        {
            phi = Newton_Raphson_Fresnel_Psi_D(
                tau->k_vals[current_point],
                tau->rho_km_vals[center],
                tau->rho_km_vals[current_point],
                tau->phi_rad_vals[current_point],
                tau->phi_rad_vals[current_point],
                tau->B_rad_vals[current_point],
                tau->EPS,
                tau->toler,
                tau->rx_km_vals[current_point],
                tau->ry_km_vals[current_point],
                tau->rz_km_vals[current_point]
            );

            x = tau->rho_km_vals[current_point] * rssringoccs_Double_Cos(phi);
            y = tau->rho_km_vals[current_point] * rssringoccs_Double_Sin(phi);
            z = tau->rz_km_vals[current_point];

            dx = x - tau->rx_km_vals[current_point];
            dy = y - tau->ry_km_vals[current_point];

            D = rssringoccs_Double_Sqrt(dx*dx + dy*dy + z*z);

            psi = -rssringoccs_Double_Fresnel_Psi(
                tau->k_vals[current_point],
                tau->rho_km_vals[center],
                tau->rho_km_vals[current_point],
                phi,
                tau->phi_rad_vals[current_point],
                tau->B_rad_vals[current_point],
                D
            );

            /*  If forward tranform is set, negate the Fresnel kernel.        */
            if (tau->use_fwd)
                psi *= -1.0;

            arg_norm = tau->window_func(window_func_x, tau->w_km_vals[center]);
            ker[i]   = rssringoccs_CDouble_Polar(arg_norm, psi);
        }
        else
            ker[i] = rssringoccs_CDouble_Zero;


        T_in[i] = tau->T_in[current_point];
    }

    fft_ker = rssringoccs_Complex_FFT(ker,  data_size, rssringoccs_False);
    if (fft_ker == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs:\n"
            "\r\trssringoccs_Diffraction_Correction_SimpleFFT\n"
            "\rrssringoccs_Complex_FFT returned NULL for fft_ker.\n"
            "\rAborting.\n\n"
        );
        return;
    }

    fft_in  = rssringoccs_Complex_FFT(T_in, data_size, rssringoccs_False);
    if (fft_in == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs:\n"
            "\r\trssringoccs_Diffraction_Correction_SimpleFFT\n"
            "\rrssringoccs_Complex_FFT returned NULL for fft_in.\n"
            "\rAborting.\n\n"
        );
        return;
    }

    for (i = 0; i < data_size; ++i)
        fft_out[i] = rssringoccs_CDouble_Multiply(fft_ker[i], fft_in[i]);

    T_out = rssringoccs_Complex_FFT(fft_out, data_size, rssringoccs_True);
    if (T_out == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs:\n"
            "\r\trssringoccs_Diffraction_Correction_SimpleFFT\n"
            "\rrssringoccs_Complex_FFT returned NULL for T_out.\n"
            "\rAborting.\n\n"
        );
        return;
    }

    for(i = 0; i < tau->n_used; ++i)
    {
        i_shift       = (nw_pts + i + shift) % (data_size);
        rcpr_F        = 1.0/tau->F_km_vals[tau->start + i];
        arg           = rssringoccs_CDouble_Rect(factor*rcpr_F, factor*rcpr_F);
        tau->T_out[i] = rssringoccs_CDouble_Multiply(arg, T_out[i_shift]);
    }

    /*  Free variables allocated by malloc.                                   */
    free(ker);
    free(fft_ker);
    free(fft_in);
    free(fft_out);
    free(T_in);
    free(T_out);
}
