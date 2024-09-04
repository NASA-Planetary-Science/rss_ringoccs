

#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_math.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_CSVData_Reverse_Geo_Variables(rssringoccs_CSVData *csv)
{
    /*  Buffer for an error message, should an error occur.                   */
    char err_mes[1024];

    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->geo)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Reverse_Geo_Variables\n\n"
            "csv->geo is NULL. Aborting\n"
        );

        return;
    }

    if (csv->geo->error_occurred)
    {
        csv->error_occurred = tmpl_True;

        /*  Keep track of error messages. Copy the previous geo message.      */
        if (csv->geo->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Reverse_Geo_Variables\n\n"
                "csv->geo has 'error_occurred' set to True.\n"
                "csv->geo set the following message:\n\n%s",
                csv->geo->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  If the function failed and no error message was set, it is likely *
         *  malloc failed. Give a generic message.                            */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_CSVData_Reverse_Geo_Variables\n\n"
                "csv->geo has 'error_occurred' set to True.\n"
            );
        }

        return;
    }

    tmpl_Double_Array_Reverse(csv->geo->rho_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->rho_dot_kms_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->D_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->rx_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->ry_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->rz_km_vals, csv->geo->n_elements);
}
