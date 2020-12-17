#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

void check_tau_data(rssringoccs_TAUObj *tau)
{
    /*  If any of these pointers are void, return 0. Else, pass.              */
    if (!tau->T_in)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has T_in set to NULL. Returning.\n\n",
            tau
        );
    else if (!tau->rho_km_vals)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has rho_km_vals set to NULL. Returning.\n\n",
            tau
        );
    else if(!tau->F_km_vals)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has F_km_vals set to NULL. Returning.\n\n",
            tau
        );
    else if (!tau->phi_rad_vals)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has phi_rad_vals set to NULL. Returning.\n\n",
            tau
        );
    else if (!tau->kd_vals)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has kd_vals set to NULL. Returning.\n\n",
            tau
        );
    else if (!tau->B_rad_vals)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has B_rad_vals set to NULL. Returning.\n\n",
            tau
        );
    else if (!tau->D_km_vals)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has D_km_vals set to NULL. Returning.\n\n",
            tau
        );
    else if (!tau->w_km_vals)
        rssringoccs_Set_Tau_Error_Message(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\tcheck_tau_data\n\n"
            "\rInput tau has w_km_vals set to NULL. Returning.\n\n",
            tau
        );
}
