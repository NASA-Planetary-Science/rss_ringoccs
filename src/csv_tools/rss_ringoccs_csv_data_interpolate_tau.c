#include <stdio.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_interpolate.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_CSVData_Interpolate_Tau(rssringoccs_CSVData *csv)
{
    char err_mes[1024];
    size_t n;

    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->tau)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Interpolate_Tau\n\n"
            "csv->tau is NULL. Aborting\n"
        );

        return;
    }

    if (csv->tau->error_occurred)
    {
        csv->error_occurred = tmpl_True;

        /*  Keep track of error messages. Copy the previous tau message.      */
        if (csv->tau->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Interpolate_Tau\n\n"
                "csv->tau has 'error_occurred' set to True.\n"
                "csv->tau set the following message:\n\n%s",
                csv->tau->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  If the function failed and no error message was set, it is likely *
         *  malloc failed. Give a generic message.                            */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Interpolate_Tau\n\n"
                "csv->tau has 'error_occurred' set to True.\n"
            );
        }

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

    if (csv->use_deprecated)
    {
        tmpl_Double_Sorted_Linear_Interp1d(
            csv->tau->rho_km_vals, csv->tau->tau_vals, csv->tau->n_elements,
            csv->rho_km_vals, csv->tau_power_vals, csv->n_elements
        );

        for (n = 0; n < csv->n_elements; ++n)
        {
            const double abs_B = tmpl_Double_Abs(csv->B_deg_vals[n]);
            const double mu = tmpl_Double_Sind(abs_B);
            const double log_power = -csv->tau_power_vals[n] / mu;
            csv->tau_power_vals[n] = tmpl_Double_Exp(log_power);
        }
    }

    else
    {
        tmpl_Double_Sorted_Linear_Interp1d(
            csv->tau->rho_km_vals, csv->tau->power_vals, csv->tau->n_elements,
            csv->rho_km_vals, csv->tau_power_vals, csv->n_elements
        );
    }
}
