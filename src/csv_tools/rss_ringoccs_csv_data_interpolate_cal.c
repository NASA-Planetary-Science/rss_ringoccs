#include <stdio.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_interpolate.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_CSVData_Interpolate_Cal(rssringoccs_CSVData *csv)
{
    char err_mes[1024];
    double *freq;
    size_t n;

    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->cal)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Interpolate_Cal\n\n"
            "csv->cal is NULL. Aborting\n"
        );

        return;
    }

    if (csv->cal->error_occurred)
    {
        csv->error_occurred = tmpl_True;

        /*  Keep track of error messages. Copy the previous cal message.      */
        if (csv->cal->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Interpolate_Cal\n\n"
                "csv->cal has 'error_occurred' set to True.\n"
                "csv->cal set the following message:\n\n%s",
                csv->cal->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  If the function failed and no error message was set, it is likely *
         *  malloc failed. Give a generic message.                            */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Interpolate_Cal\n\n"
                "csv->cal has 'error_occurred' set to True.\n"
            );
        }

        return;
    }

    /*  Steal the pointer to avoid a call to malloc. It is destroyed in the   *
     *  end anyways, so no harm done.                                         */
    freq = csv->cal->f_sky_pred_vals;

    /*  The sky frequency is the difference of the predicted and residual     *
     *  frequencies. Loop through the arrays and compute the difference.      */
    for (n = 0; n < csv->cal->n_elements; ++n)
        freq[n] -= csv->cal->f_sky_resid_fit_vals[n];

    tmpl_Double_Sorted_Interp1d(
        csv->cal->t_oet_spm_vals, freq, csv->cal->n_elements,
        csv->t_oet_spm_vals, csv->f_sky_hz_vals, csv->n_elements
    );
}
