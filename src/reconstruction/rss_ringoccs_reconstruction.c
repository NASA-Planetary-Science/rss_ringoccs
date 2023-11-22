
#include <stdlib.h>
#include <stdio.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau)
{
    tmpl_ComplexDouble *temp_T_in;
    tmpl_Bool temp_fwd;
    size_t n, temp_start, temp_n_used, nw_pts;
    double w_left, w_right, w_max;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    rssringoccs_Tau_Check_Keywords(tau);
    rssringoccs_Tau_Check_Occ_Type(tau);
    rssringoccs_Tau_Get_Window_Width(tau);
    rssringoccs_Tau_Check_Data_Range(tau);

    tau->T_out = calloc(tau->arr_size, sizeof(*tau->T_out));
    rssringoccs_Tau_Check_Data(tau);

    temp_fwd = tau->use_fwd;
    tau->use_fwd = tmpl_False;

    if      (tau->psinum == rssringoccs_DR_Fresnel)
        rssringoccs_Diffraction_Correction_Fresnel(tau);
    else if (tau->psinum == rssringoccs_DR_Legendre)
        rssringoccs_Diffraction_Correction_Legendre(tau);
    else if (tau->psinum == rssringoccs_DR_NewtonSimpleFFT)
        rssringoccs_Diffraction_Correction_SimpleFFT(tau);
    else
        rssringoccs_Diffraction_Correction_Newton(tau);

    tau->use_fwd = temp_fwd;

    if (tau->use_fwd)
    {
        temp_T_in  = tau->T_in;
        tau->T_in  = tau->T_out;
        tau->T_out = calloc(tau->arr_size, sizeof(*tau->T_out));

        /*  If forward tranform is set, negate the k_vals variable. This has  *
         *  the equivalent effect of computing the forward calculation later. */
        for (n = 0; n <= tau->n_used; ++n)
            tau->k_vals[tau->start + n] *= -1.0;

        w_left  = tau->w_km_vals[tau->start];
        w_right = tau->w_km_vals[tau->start + tau->n_used];

        if (w_left < w_right)
            w_max = w_right;
        else
            w_max = w_left;

        nw_pts = (size_t) (w_max / (tau->dx_km * 2.0));
        temp_start = tau->start;
        temp_n_used = tau->n_used;

        if (tau->n_used <= 2*nw_pts)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_strdup(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Reconstruction\n\n"
                "\rNot enough data available to perform the forward model.\n"
                "\rReturning with T_fwd pointer set to an array of zeroes.\n"
            );
        }
        else
        {
            tau->start = tau->start + nw_pts;
            tau->n_used = tau->n_used - 2*nw_pts;
            if      (tau->psinum == rssringoccs_DR_Fresnel)
                rssringoccs_Diffraction_Correction_Fresnel(tau);
            else if (tau->psinum == rssringoccs_DR_Legendre)
                rssringoccs_Diffraction_Correction_Legendre(tau);
            else if (tau->psinum == rssringoccs_DR_NewtonSimpleFFT)
                rssringoccs_Diffraction_Correction_SimpleFFT(tau);
            else
                rssringoccs_Diffraction_Correction_Newton(tau);

            tau->start = temp_start;
            tau->n_used = temp_n_used;
        }

        /*  Unnegate k_vals to its original values.                           */
        for (n = 0; n <= tau->n_used; ++n)
            tau->k_vals[tau->start + n] *= -1.0;

        tau->T_fwd = tau->T_out;
        tau->T_out = tau->T_in;
        tau->T_in  = temp_T_in;

    }

    rssringoccs_Tau_Finish(tau);
    return;
}
