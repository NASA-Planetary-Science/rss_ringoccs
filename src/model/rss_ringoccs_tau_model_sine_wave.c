#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_complex.h>
#include <libtmpl/include/tmpl_math.h>
#include <rss_ringoccs/include/types/rss_ringoccs_tauobj.h>
#include <rss_ringoccs/include/rss_ringoccs_model.h>
#include <stddef.h>

void
rssringoccs_Tau_Model_Sine_Wave(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    const rssringoccs_ModelParameters * TMPL_RESTRICT const parameters
)
{
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
            "\r\trssringoccs_Tau_Model_Sine_Wave\n\n"
            "\rparameters pointer is NULL.\n\n";

        return;
    }

    if (parameters->model != rssringoccs_Model_SineWave)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Model_Sine_Wave\n\n"
            "\rModel type is not set to rssringoccs_Model_SineWave.\n\n";

        return;
    }

    for (n = 0; n < tau->dlp->arr_size; ++n)
    {
        const double x = tau->dlp->rho_km_vals[n];
        const double t = parameters->geometry.wave.frequency * x;
        const double y = tmpl_Double_Sin(t);
        const double p = 0.5 * parameters->peak_opacity * (1.0 + y);
        tau->T_in[n] = tmpl_CDouble_Rect(p, 0.0);
    }
}
