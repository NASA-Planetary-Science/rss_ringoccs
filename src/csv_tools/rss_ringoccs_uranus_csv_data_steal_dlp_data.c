
/*  Functions for reading and writing files.                                  */
#include <stdio.h>

/*  malloc provided here.                                                     */
#include <stdlib.h>

/*  libtmpl provides Booleans, string duplicate, and math tools.              */
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_math.h>

/*  Prototype for the function and typedefs for structs.                      */
#include <rss_ringoccs/include/rss_ringoccs_csv_tools.h>

void rssringoccs_UranusCSVData_Steal_DLP_Data(rssringoccs_UranusCSVData *csv)
{
    /*  Buffer for an error message, should an error occur.                   */
    char err_mes[1024];

    if (!csv)
        return;

    if (csv->error_occurred)
        return;

    if (!csv->dlp)
    {
        csv->error_occurred = tmpl_True;
        csv->error_message = tmpl_String_Duplicate(
            "\nError Encountered: rss_ringoccs\n"
            "\trssringoccs_UranusCSVData_Steal_DLP_Data\n\n"
            "csv->dlp is NULL. Aborting\n"
        );

        return;
    }

    if (csv->dlp->error_occurred)
    {
        csv->error_occurred = tmpl_True;

        /*  Keep track of error messages. Copy the previous geo message.      */
        if (csv->dlp->error_message)
        {
            sprintf(
                err_mes,
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Steal_DLP_Data\n\n"
                "csv->dlp has 'error_occurred' set to True.\n"
                "csv->dlp set the following message:\n\n%s",
                csv->dlp->error_message
            );

            csv->error_message = tmpl_String_Duplicate(err_mes);
        }

        /*  If the function failed and no error message was set, it is likely *
         *  malloc failed. Give a generic message.                            */
        else
        {
            csv->error_message = tmpl_String_Duplicate(
                "\nError Encountered: rss_ringoccs\n"
                "\trssringoccs_UranusCSVData_Steal_DLP_Data\n\n"
                "csv->dlp has 'error_occurred' set to True.\n"
            );
        }

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
    csv->phi_deg_vals = csv->dlp->phi_deg_vals;
    csv->B_deg_vals = csv->dlp->B_deg_vals;
    csv->phi_rl_deg_vals = csv->dlp->phi_rl_deg_vals;

    /*  The new TAB files have normalized power included.                     */
    csv->p_norm_vals = csv->dlp->p_norm_vals;

    /*  The Uranus DLP files have more data than the Saturn DLPs. This makes  *
     *  the CAL files not needed. Steal this data too.                        */
    csv->rho_dot_kms_vals = csv->dlp->rho_dot_kms_vals;
    csv->D_km_vals = csv->dlp->D_km_vals;
    csv->f_sky_hz_vals = csv->dlp->f_sky_hz_vals;
}
