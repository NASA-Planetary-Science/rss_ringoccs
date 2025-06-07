
#include <stddef.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/tmpl_math.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_model.h>

void
rssringoccs_Tau_Model_Left_Straightedge(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const rssringoccs_ModelParameters * TMPL_RESTRICT const parameters
)
{
    tmpl_ComplexDouble transmittance;
    size_t n, nw_pts;
    double w_max;
    tmpl_Bool use_fwd;

    if (!tau)
        return;

    if (tau->error_occurred)
        return;

    if (!parameters)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Model_Left_Straightedge\n\n"
            "\rparameters pointer is NULL.\n\n";

        return;
    }

    if (parameters->model != rssringoccs_Model_LeftEdge)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Model_Left_Straightedge\n\n"
            "\rModel type is not set to rssringoccs_Model_LeftEdge.\n\n";

        return;
    }

    tau->T_out = TMPL_MALLOC(tmpl_ComplexDouble, tau->arr_size);

    rssringoccs_Tau_Check_Core_Data(tau);

    if (tau->error_occurred)
        return;

    transmittance = tmpl_CDouble_Rect(parameters->peak_opacity, 0.0);

    if (parameters->geometry.edge.center <= tau->rho_km_vals[0])
    {
        for (n = 0; n < tau->arr_size; ++n)
            tau->T_in[n] = tmpl_CDouble_Zero;

        return;
    }

    if (tau->rho_km_vals[tau->arr_size - 1] <= parameters->geometry.edge.center)
    {
        for (n = 0; n < tau->arr_size; ++n)
            tau->T_in[n] = transmittance;

        return;
    }

    n = 0;

    while (tau->rho_km_vals[n] < parameters->geometry.edge.center)
    {
        tau->T_in[n] = transmittance;
        ++n;
    }

    for (; n < tau->arr_size; ++n)
        tau->T_in[n] = tmpl_CDouble_Zero;

    use_fwd = tau->use_fwd;
    tau->use_fwd = tmpl_True;
    rssringoccs_Diffraction_Correction(tau);

    for (n = 0; n < tau->arr_size; ++n)
    {
        tau->T_in[n] = tau->T_out[n];
        tau->T_out[n] = tmpl_CDouble_Zero;
    }

    w_max = tmpl_Double_Array_Max(tau->w_km_vals, tau->arr_size);

    nw_pts = TMPL_CAST(w_max / (tau->dx_km * 2.0), size_t);

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
    tau->use_fwd = tmpl_False;

    rssringoccs_Diffraction_Correction(tau);

    tau->use_fwd = use_fwd;

    if (tau->use_fwd)
    {
        tmpl_ComplexDouble * const temp_T_in = tau->T_in;
        tau->T_in = tau->T_out;
        tau->T_out = calloc(tau->arr_size, sizeof(*tau->T_out));

        nw_pts = (size_t) (w_max / (tau->dx_km * 2.0));

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

        tau->T_fwd = tau->T_out;
        tau->T_out = tau->T_in;
        tau->T_in  = temp_T_in;
    }
}
