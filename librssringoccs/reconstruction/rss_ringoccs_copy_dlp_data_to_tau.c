#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_calibration.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stdlib.h>

void
rssringoccs_Copy_DLP_Data_To_Tau(rssringoccs_DLPObj *dlp,
                                 rssringoccs_TAUObj *tau)
{
    unsigned long n;
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (dlp == NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "Input dlp pointer is NULL. Returning.\n"
        );
        return;
    }

    if (dlp->error_occurred)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Copy_DLP_Data_To_Tau\n\n"
            "Input dlp pointer has the error_occurred member set to True.\n"
        );
        return;
    }

    tau->arr_size = dlp->arr_size;
    tau->rho_km_vals = malloc(sizeof(*tau->rho_km_vals) * tau->arr_size);
    tau->phi_rad_vals = malloc(sizeof(*tau->phi_rad_vals) * tau->arr_size);
    tau->f_sky_hz_vals = malloc(sizeof(*tau->f_sky_hz_vals) * tau->arr_size);
    tau->rho_dot_kms_vals
        = malloc(sizeof(*tau->rho_dot_kms_vals) * tau->arr_size);

    tau->raw_tau_threshold_vals
        = malloc(sizeof(*tau->raw_tau_threshold_vals) * tau->arr_size);

    tau->B_rad_vals = malloc(sizeof(*tau->B_rad_vals) * tau->arr_size);
    tau->D_km_vals = malloc(sizeof(*tau->D_km_vals) * tau->arr_size);
    tau->t_oet_spm_vals = malloc(sizeof(*tau->t_oet_spm_vals) * tau->arr_size);
    tau->t_ret_spm_vals = malloc(sizeof(*tau->t_ret_spm_vals) * tau->arr_size);
    tau->t_set_spm_vals = malloc(sizeof(*tau->t_set_spm_vals) * tau->arr_size);

    tau->rho_corr_pole_km_vals
        = malloc(sizeof(*tau->rho_corr_pole_km_vals) * tau->arr_size);

    tau->rho_corr_timing_km_vals
        = malloc(sizeof(*tau->rho_corr_timing_km_vals) * tau->arr_size);

    tau->phi_rl_rad_vals
        = malloc(sizeof(*tau->phi_rl_rad_vals) * tau->arr_size);

    tau->p_norm_vals = malloc(sizeof(*tau->p_norm_vals) * tau->arr_size);
    tau->phase_rad_vals = malloc(sizeof(*tau->phase_rad_vals) * tau->arr_size);

    for (n=0; n<dlp->arr_size; ++n)
    {
        tau->rho_km_vals[n] = dlp->rho_km_vals[n];
        tau->phi_rad_vals[n] = dlp->phi_rad_vals[n];
        tau->B_rad_vals[n] = dlp->B_rad_vals[n];
        tau->D_km_vals[n] = dlp->D_km_vals[n];
        tau->f_sky_hz_vals[n] = dlp->f_sky_hz_vals[n];
        tau->rho_dot_kms_vals[n] = dlp->rho_dot_kms_vals[n];
        tau->t_oet_spm_vals[n] = dlp->t_oet_spm_vals[n];
        tau->t_ret_spm_vals[n] = dlp->t_ret_spm_vals[n];
        tau->t_set_spm_vals[n] = dlp->t_set_spm_vals[n];
        tau->rho_corr_pole_km_vals[n] = dlp->rho_corr_pole_km_vals[n];
        tau->rho_corr_timing_km_vals[n] = dlp->rho_corr_timing_km_vals[n];
        tau->phi_rl_rad_vals[n] = dlp->phi_rl_rad_vals[n];
        tau->p_norm_vals[n] = dlp->p_norm_vals[n];

        /*  The phase needs to be negated due to mathematical conventions.    */
        tau->phase_rad_vals[n] = -dlp->phase_rad_vals[n];
        tau->raw_tau_threshold_vals[n] = dlp->raw_tau_threshold_vals[n];
    }

}
