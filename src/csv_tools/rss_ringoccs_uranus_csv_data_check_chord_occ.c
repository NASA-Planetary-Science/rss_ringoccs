#include <stdio.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_UranusCSVData_Check_Chord_Occ(rssringoccs_UranusCSVData *csv)
{
    char err_mes[1024];
    double min_dr_dt, max_dr_dt;
    size_t n;

    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->geo)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Check_Chord_Occ\n\n"
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
                "\trssringoccs_UranusCSVData_Check_Chord_Occ\n\n"
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
                "\trssringoccs_UranusCSVData_Check_Chord_Occ\n\n"
                "csv->geo has 'error_occurred' set to True.\n"
            );
        }

        return;
    }

    /*  Since "dt" is about on the order of 1 second, and since the width of  *
     *  the universe is ~10^23 km, the ratio should not be too large. We can  *
     *  safely initialize the min and max to +/- 10^36, which portably fits   *
     *  in both single and double precision numbers, and avoid an "infinity"  *
     *  call.                                                                 */
    min_dr_dt = +1.0E+36;
    max_dr_dt = -1.0E+36;

    /*  Find the minimum and maximum of drho/dt.                              */
    for (n = 0; n < csv->n_elements - 1; ++n)
    {
        /*  Compute the very simple numerical derivative using the samples.   */
        const double dr = csv->rho_km_vals[n+1U] - csv->rho_km_vals[n];
        const double dt = csv->t_set_spm_vals[n+1U] - csv->t_set_spm_vals[n];
        const double deriv = dr / dt;

        /*  Search for the min and max of the derivative.                     */
        if (deriv < min_dr_dt)
            min_dr_dt = deriv;

        if (max_dr_dt < deriv)
            max_dr_dt = deriv;
    }

    /*  Check for errors.                                                     */
    if ((min_dr_dt < 0.0) && (max_dr_dt > 0.0))
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Check_Chord_Occ\n\n"
            "\rdrho/dt has positive and negative values. Check your DLP file.\n"
            "\rIt is likely a chord occultation and needs to be split into\n"
            "\ringress and egress portions.\n"
        );

        return;
    }

    if ((min_dr_dt == 0.0) || (max_dr_dt == 0.0))
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "Error Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Check_Chord_Occ\n\n"
            "\rdrho/dt has zero-valued elements. Check your DLP file.\n"
            "\rIt is likely a chord occultation and needs to be split into\n"
            "\ringress and egress portions.\n"
        );

        return;
    }

    if (max_dr_dt < 0.0)
    {
        rssringoccs_UranusCSVData_Reverse_Geo_Variables(csv);

        while (csv->geo->rho_dot_kms_vals[csv->geo_increment] >= 0.0)
        {
            if (csv->geo_increment >= csv->geo->n_elements)
                break;

            ++csv->geo_increment;
        }
    }

    else if (min_dr_dt > 0.0)
    {
        while (csv->geo->rho_dot_kms_vals[csv->geo_increment] <= 0.0)
        {
            if (csv->geo_increment >= csv->geo->n_elements)
                break;

            ++csv->geo_increment;
        }
    }
}
