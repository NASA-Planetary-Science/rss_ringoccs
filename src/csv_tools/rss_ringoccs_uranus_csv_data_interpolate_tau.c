#include <stdio.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_interpolate.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_UranusCSVData_Interpolate_Tau(rssringoccs_UranusCSVData *csv)
{
    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->tau)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Interpolate_Tau\n\n"
            "csv->tau is NULL.\n\n";

        return;
    }

    if (csv->tau->error_occurred)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Interpolate_Tau\n\n"
            "csv->tau has error_occurred set to True.\n\n";

        return;
    }

    tmpl_Double_Sorted_Linear_Interp1d(
        csv->tau->rho_km_vals, csv->tau->phase_deg_vals, csv->tau->n_elements,
        csv->rho_km_vals, csv->tau_phase_deg_vals, csv->n_elements
    );

    tmpl_Double_Sorted_Linear_Interp1d(
        csv->tau->rho_km_vals, csv->tau->tau_vals, csv->tau->n_elements,
        csv->rho_km_vals, csv->tau_vals, csv->n_elements
    );

    tmpl_Double_Sorted_Linear_Interp1d(
        csv->tau->rho_km_vals, csv->tau->power_vals, csv->tau->n_elements,
        csv->rho_km_vals, csv->tau_power_vals, csv->n_elements
    );
}
