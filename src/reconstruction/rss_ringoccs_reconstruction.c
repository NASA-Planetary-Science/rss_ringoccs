#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/compat/tmpl_calloc.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stdlib.h>

#if 1

void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau)
{
    tmpl_Bool temp_fwd;
    size_t temp_start, temp_n_used, nw_pts;
    double w_left, w_right, w_max;

    if (!tau)
        return;

    rssringoccs_DLP_Check_Core_Data(tau->dlp);
    rssringoccs_DLP_Check_Geometry(tau->dlp);
    rssringoccs_DLP_Check_Occ_Type(tau->dlp);

    if (tau->dlp->error_occurred)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tau->dlp->error_message;
    }

    rssringoccs_Tau_Check_Keywords(tau);
    rssringoccs_Tau_Get_Window_Width(tau);

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Core_Data(tau);

    /*  Check to ensure you have enough data to process.                      */
    rssringoccs_Tau_Check_Data_Range(tau);

    /*  The previous functions set the error_occurred Boolean on failure.     */
    if (tau->error_occurred)
        return;

    tau->T_out = TMPL_MALLOC(tmpl_ComplexDouble, tau->dlp->arr_size);

    temp_fwd = tau->use_fwd;
    tau->use_fwd = tmpl_False;

    rssringoccs_Diffraction_Correction(tau);

    tau->use_fwd = temp_fwd;

    if (tau->use_fwd)
    {
        tmpl_ComplexDouble * const temp_T_in = tau->T_in;
        tau->T_in = tau->T_out;
        tau->T_out = TMPL_CALLOC(tmpl_ComplexDouble, tau->dlp->arr_size);

        w_left  = tau->w_km_vals[tau->start];
        w_right = tau->w_km_vals[tau->start + tau->n_used];

        if (w_left < w_right)
            w_max = w_right;
        else
            w_max = w_left;

        nw_pts = (size_t) (w_max / (2.0 * tau->dlp->dx_km));
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

    return;
}

#else

#include <rss_ringoccs/include/rss_ringoccs_model.h>

void rssringoccs_Reconstruction(rssringoccs_TAUObj *tau)
{
    rssringoccs_ModelParameters parameters;

    if (!tau)
        return;

    rssringoccs_Tau_Check_Keywords(tau);
    rssringoccs_Tau_Check_Occ_Type(tau);
    rssringoccs_Tau_Get_Window_Width(tau);
    rssringoccs_Tau_Check_Data_Range(tau);

    parameters.model = rssringoccs_Model_LeftEdge;
    parameters.peak_opacity = 1.0;
    parameters.geometry.edge.center = 87500.0;

    rssringoccs_Tau_Model_Left_Straightedge(tau, &parameters);

    return;
}

#endif
