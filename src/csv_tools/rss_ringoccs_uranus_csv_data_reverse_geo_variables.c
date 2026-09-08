#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void
rssringoccs_UranusCSVData_Reverse_Geo_Variables(rssringoccs_UranusCSVData *csv)
{
    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->geo)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Reverse_Geo_Variables\n\n"
            "csv->geo is NULL.\n\n";

        return;
    }

    if (csv->geo->error_occurred)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message =
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Reverse_Geo_Variables\n\n"
            "csv->geo has error_occurred set to True.\n\n";

        return;
    }

    tmpl_Double_Array_Reverse(csv->geo->rho_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->rho_dot_kms_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->D_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->rx_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->ry_km_vals, csv->geo->n_elements);
    tmpl_Double_Array_Reverse(csv->geo->rz_km_vals, csv->geo->n_elements);
}
