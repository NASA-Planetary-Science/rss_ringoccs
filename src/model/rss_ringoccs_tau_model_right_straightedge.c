
#include <stddef.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_model.h>

void
rssringoccs_Tau_Model_Right_Straightedge(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const rssringoccs_ModelParameters * TMPL_RESTRICT const parameters
)
{
    tmpl_ComplexDouble transmittance;
    size_t n;

    if (!tau)
        return;

    if (tau->error_occurred)
        return;

    if (!parameters)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Model_Right_Straightedge\n\n"
            "\rparameters pointer is NULL.\n\n";

        return;
    }

    if (parameters->model != rssringoccs_Model_RightEdge)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Model_Right_Straightedge\n\n"
            "\rModel type is not set to rssringoccs_Model_RightEdge.\n\n";

        return;
    }

    transmittance = tmpl_CDouble_Rect(parameters->peak_opacity, 0.0);

    if (parameters->geometry.edge.center <= tau->rho_km_vals[0])
    {
        for (n = 0; n < tau->arr_size; ++n)
            tau->T_in[n] = transmittance;

        return;
    }

    if (tau->rho_km_vals[tau->arr_size - 1] <= parameters->geometry.edge.center)
    {
        for (n = 0; n < tau->arr_size; ++n)
            tau->T_in[n] = tmpl_CDouble_Zero;

        return;
    }

    n = 0;

    while (tau->rho_km_vals[n] < parameters->geometry.edge.center)
    {
        tau->T_in[n] = tmpl_CDouble_Zero;
        ++n;
    }

    for (; n < tau->arr_size; ++n)
        tau->T_in[n] = transmittance;
}
