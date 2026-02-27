#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>
#include <rss_ringoccs/include/rss_ringoccs_model.h>
#include <stddef.h>

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

    for (n = 0; n < tau->dlp->arr_size; ++n)
    {
        if (tau->dlp->rho_km_vals[n] < parameters->geometry.edge.center)
            tau->T_in[n] = tmpl_CDouble_Zero;
        else
            tau->T_in[n] = transmittance;
    }
}
