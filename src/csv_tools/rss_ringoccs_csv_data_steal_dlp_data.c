/*  malloc provided here.                                                     */
#include <stdlib.h>

/*  libtmpl provides Booleans, string duplicate, and math tools.              */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_CSVData_Steal_DLP_Data(rssringoccs_CSVData *csv)
{
    size_t n;
    const size_t zero = (size_t)0;

    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->dlp)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Steal_DLP_Data\n\n"
            "csv->dlp is NULL.\n\n";

        return;
    }

    if (csv->dlp->error_occurred)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Steal_DLP_Data\n\n"
            "csv->dlp has error_occurred set to True.\n\n";

        return;
    }

    /*  The references to these variables can be stolen to avoid copying.     */
    csv->rho_km_vals = csv->dlp->rho_km_vals;
    csv->raw_tau_vals = csv->dlp->raw_tau_vals;
    csv->t_ret_spm_vals = csv->dlp->t_ret_spm_vals;
    csv->t_set_spm_vals = csv->dlp->t_set_spm_vals;
    csv->t_oet_spm_vals = csv->dlp->t_oet_spm_vals;
    csv->rho_corr_pole_km_vals = csv->dlp->rho_corr_pole_km_vals;
    csv->rho_corr_timing_km_vals = csv->dlp->rho_corr_timing_km_vals;
    csv->raw_tau_threshold_vals = csv->dlp->raw_tau_threshold_vals;

    /*  Various angles from the dlp file. Steal the reference.                */
    csv->phase_deg_vals = csv->dlp->phase_deg_vals;
    csv->phi_deg_vals = csv->dlp->phi_ora_deg_vals;
    csv->B_deg_vals = csv->dlp->B_deg_vals;
    csv->phi_rl_deg_vals = csv->dlp->phi_rl_deg_vals;

    /*  The new TAB files have normalized power included.                     */
    if (!csv->use_deprecated)
        csv->p_norm_vals = csv->dlp->p_norm_vals;

    /*  Otherwise we need to compute the power from the optical depth.        */
    else
    {
        /*  Need to allocate memory for p_norm_vals.                          */
        csv->p_norm_vals = malloc(sizeof(*csv->p_norm_vals) * csv->n_elements);

        if (!csv->p_norm_vals)
        {
            csv->error_occurred = tmpl_True;
            csv->error_message =
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Steal_DLP_Data\n\n"
                "malloc returned NULL for p_norm_vals.\n\n";

            return;
        }

        /*  Loop through the elements of the DLP tau variable and compute the *
         *  normalize diffracted power.                                       */
        for (n = zero; n < csv->n_elements; ++n)
        {
            /*  tmpl_Double_Sind computes sine of an argument in degrees.     */
            const double B = csv->dlp->B_deg_vals[n];
            const double mu = tmpl_Double_Sind(tmpl_Double_Abs(B));
            const double log_power = -csv->raw_tau_vals[n] / mu;

            /*  Normalize diffracted power can be computed from optical depth.*/
            csv->p_norm_vals[n] = tmpl_Double_Exp(log_power);
        }
    }
}
