#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_interpolate.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_CSVData_Interpolate_Geo(rssringoccs_CSVData *csv)
{
    char err_mes[1024];
    double *rho, *rho_dot, *D, *rx, *ry, *rz;
    size_t length, decrease;

    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->geo)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Interpolate_Geo\n\n"
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
                "\trssringoccs_CSVData_Interpolate_Geo\n\n"
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
                "\trssringoccs_CSVData_Interpolate_Geo\n\n"
                "csv->geo has 'error_occurred' set to True.\n"
            );
        }

        return;
    }

    decrease = csv->geo_increment + csv->geo_decrement;

    if (csv->geo->n_elements < decrease)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_CSVData_Interpolate_Geo\n\n"
            "csv->geo->n_elements less than sum of increment and decrement.\n"
        );

        return;
    }

    length = csv->geo->n_elements - decrease;

    rho = csv->geo->rho_km_vals + csv->geo_increment;
    rho_dot = csv->geo->rho_dot_kms_vals + csv->geo_increment;
    D = csv->geo->D_km_vals + csv->geo_increment;
    rx = csv->geo->rx_km_vals + csv->geo_increment;
    ry = csv->geo->ry_km_vals + csv->geo_increment;
    rz = csv->geo->rz_km_vals + csv->geo_increment;

    tmpl_Double_Sorted_Interp1d(
        rho, D, length,
        csv->rho_km_vals, csv->D_km_vals, csv->n_elements
    );

    tmpl_Double_Sorted_Interp1d(
        rho, rho_dot, length,
        csv->rho_km_vals, csv->rho_dot_kms_vals, csv->n_elements
    );

    tmpl_Double_Sorted_Interp1d(
        rho, rx, length,
        csv->rho_km_vals, csv->rx_km_vals, csv->n_elements
    );

    tmpl_Double_Sorted_Interp1d(
        rho, ry, length,
        csv->rho_km_vals, csv->ry_km_vals, csv->n_elements
    );

    tmpl_Double_Sorted_Interp1d(
        rho, rz, length,
        csv->rho_km_vals, csv->rz_km_vals, csv->n_elements
    );
}
