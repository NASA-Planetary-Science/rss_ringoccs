#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>

void rssringoccs_Check_Tau_Data(rssringoccs_TAUObj *tau)
{
    if (tau == NULL)
        return;

    /*  If any of these pointers are void, return 0. Else, pass.              */
    if (tau->T_in == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Check_Tau_Data\n\n"
            "\rInput tau has T_in set to NULL. Returning.\n\n"
        );
    }
    else if (tau->rho_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has rho_km_vals set to NULL. Returning.\n\n"
        );
    }
    else if (tau->F_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has F_km_vals set to NULL. Returning.\n\n"
        );
    }
    else if (tau->phi_rad_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has phi_rad_vals set to NULL. Returning.\n\n"
        );
    }
    else if (tau->kd_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has kd_vals set to NULL. Returning.\n\n"
        );
    }
    else if (tau->B_rad_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has B_rad_vals set to NULL. Returning.\n\n"
        );
    }
    else if (tau->D_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has D_km_vals set to NULL. Returning.\n\n"
        );
    }
    else if (tau->w_km_vals == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has w_km_vals set to NULL. Returning.\n\n"
        );
    }
}
