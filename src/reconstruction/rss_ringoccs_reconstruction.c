#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stdlib.h>

void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau)
{
    tmpl_Bool temp_fwd;
    size_t temp_start, temp_n_used, nw_pts;
    double w_left, w_right, w_max;

    if (!tau)
        return;

    rssringoccs_Tau_Check_Keywords(tau);
    rssringoccs_Tau_Check_Occ_Type(tau);
    rssringoccs_Tau_Get_Window_Width(tau);
    rssringoccs_Tau_Check_Data_Range(tau);

    if (tau->error_occurred)
        return;

    tau->T_out = TMPL_MALLOC(tmpl_ComplexDouble, tau->arr_size);

    rssringoccs_Tau_Check_Core_Data(tau);

    if (tau->error_occurred)
        return;

    temp_fwd = tau->use_fwd;
    tau->use_fwd = tmpl_False;

    rssringoccs_Diffraction_Correction(tau);

    tau->use_fwd = temp_fwd;

    if (tau->use_fwd)
    {
        tmpl_ComplexDouble * const temp_T_in = tau->T_in;
        tau->T_in = tau->T_out;
        tau->T_out = calloc(tau->arr_size, sizeof(*tau->T_out));

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
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Reconstruction\n\n"
                "\rNot enough data available to perform the forward model.\n"
                "\rReturning with T_fwd pointer set to an array of zeroes.\n";

            return;
        }

        tau->start = tau->start + nw_pts;
        tau->n_used = tau->n_used - 2*nw_pts;

        rssringoccs_Diffraction_Correction(tau);

        tau->start = temp_start;
        tau->n_used = temp_n_used;
        tau->T_fwd = tau->T_out;
        tau->T_out = tau->T_in;
        tau->T_in  = temp_T_in;
    }

    rssringoccs_Tau_Finish(tau);
    return;
}
