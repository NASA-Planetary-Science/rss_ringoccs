
#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

static void __resize_array(double **ptr, unsigned long start, unsigned long len)
{
    double *temp, *data;
    unsigned long n;
    temp = malloc(sizeof(*temp) * len);

    data = *ptr;

    for (n = 0; n < len; ++n)
        temp[n] = data[start + n];

    free(data);
    *ptr = temp;
}

void rssringoccs_Tau_Finish(rssringoccs_TAUObj* tau)
{
    double mu, factor;
    unsigned long n;
    unsigned long len;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    len = tau->n_used;

    tau->power_vals         = malloc(sizeof(*tau->power_vals)         * len);
    tau->phase_vals         = malloc(sizeof(*tau->phase_vals)         * len);
    tau->tau_vals           = malloc(sizeof(*tau->tau_vals)           * len);
    tau->tau_threshold_vals = malloc(sizeof(*tau->tau_threshold_vals) * len);

    __resize_array(&tau->rho_km_vals, tau->start, len);
    __resize_array(&tau->F_km_vals, tau->start, len);
    __resize_array(&tau->phi_rad_vals, tau->start, len);
    __resize_array(&tau->k_vals, tau->start, len);
    __resize_array(&tau->f_sky_hz_vals, tau->start, len);
    __resize_array(&tau->rho_dot_kms_vals, tau->start, len);
    __resize_array(&tau->raw_tau_threshold_vals, tau->start, len);
    __resize_array(&tau->B_rad_vals, tau->start, len);
    __resize_array(&tau->D_km_vals, tau->start, len);
    __resize_array(&tau->w_km_vals, tau->start, len);
    __resize_array(&tau->t_oet_spm_vals, tau->start, len);
    __resize_array(&tau->t_ret_spm_vals, tau->start, len);
    __resize_array(&tau->t_set_spm_vals, tau->start, len);
    __resize_array(&tau->rho_corr_pole_km_vals, tau->start, len);
    __resize_array(&tau->rho_corr_timing_km_vals, tau->start, len);
    __resize_array(&tau->phi_rl_rad_vals, tau->start, len);
    __resize_array(&tau->p_norm_vals, tau->start, len);
    __resize_array(&tau->phase_rad_vals, tau->start, len);

    factor = rssringoccs_Double_Log(tau->dx_km / tau->res);

    for (n = 0; n < len; ++n)
    {
        tau->power_vals[n] = rssringoccs_CDouble_Abs_Squared(tau->T_out[n]);
        tau->phase_vals[n] = rssringoccs_CDouble_Argument(tau->T_out[n]);

        mu = rssringoccs_Double_Sin(rssringoccs_Double_Abs(tau->B_rad_vals[n]));
        tau->tau_vals[n] = -mu*rssringoccs_Double_Log(tau->power_vals[n]);
        tau->tau_threshold_vals[n] = tau->raw_tau_threshold_vals[n] - factor*mu;
    }
    tau->arr_size = len;
}
