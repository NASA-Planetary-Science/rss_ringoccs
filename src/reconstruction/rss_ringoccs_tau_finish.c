
#include <stdlib.h>
#include <math.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

static void resize_array(double **ptr, size_t start, size_t len)
{
    double *temp, *data;
    size_t n;
    temp = malloc(sizeof(*temp) * len);

    data = *ptr;

    for (n = 0; n < len; ++n)
        temp[n] = data[start + n];

    free(data);
    *ptr = temp;
}

static void resize_carray(tmpl_ComplexDouble **ptr, size_t start, size_t len)
{
    tmpl_ComplexDouble *temp, *data;
    size_t n;
    temp = malloc(sizeof(*temp) * len);
    data = *ptr;

    for (n = 0; n < len; ++n)
        temp[n] = data[start + n];

    free(data);
    *ptr = temp;
}

#if 1

void rssringoccs_Tau_Finish(rssringoccs_TAUObj* tau)
{
    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    resize_carray(&tau->T_in, tau->start, tau->n_used);
    resize_carray(&tau->T_out, tau->start, tau->n_used);
    resize_array(&tau->rho_km_vals, tau->start, tau->n_used);
    resize_array(&tau->F_km_vals, tau->start, tau->n_used);
    resize_array(&tau->phi_deg_vals, tau->start, tau->n_used);
    resize_array(&tau->k_vals, tau->start, tau->n_used);
    resize_array(&tau->rho_dot_kms_vals, tau->start, tau->n_used);
    resize_array(&tau->B_deg_vals, tau->start, tau->n_used);
    resize_array(&tau->D_km_vals, tau->start, tau->n_used);
    resize_array(&tau->w_km_vals, tau->start, tau->n_used);
    resize_array(&tau->t_oet_spm_vals, tau->start, tau->n_used);
    resize_array(&tau->t_ret_spm_vals, tau->start, tau->n_used);
    resize_array(&tau->t_set_spm_vals, tau->start, tau->n_used);
    resize_array(&tau->rho_corr_pole_km_vals, tau->start, tau->n_used);
    resize_array(&tau->rho_corr_timing_km_vals, tau->start, tau->n_used);
    resize_array(&tau->phi_rl_deg_vals, tau->start, tau->n_used);
    tau->arr_size = tau->n_used;

    if (tau->use_fwd)
        resize_carray(&tau->T_fwd, tau->start, tau->n_used);
}

#else

void rssringoccs_Tau_Finish(rssringoccs_TAUObj* tau)
{
    double mu, factor;
    size_t n, len;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    len = tau->n_used;

    tau->power_vals         = malloc(sizeof(*tau->power_vals)         * len);
    tau->phase_vals         = malloc(sizeof(*tau->phase_vals)         * len);
    tau->tau_vals           = malloc(sizeof(*tau->tau_vals)           * len);
    tau->tau_threshold_vals = malloc(sizeof(*tau->tau_threshold_vals) * len);

    if (tau->use_fwd)
    {
        tau->p_norm_fwd_vals = malloc(sizeof(*tau->p_norm_fwd_vals) * len);
        tau->phase_fwd_vals  = malloc(sizeof(*tau->phase_fwd_vals)  * len);
        tau->tau_fwd_vals    = malloc(sizeof(*tau->tau_fwd_vals)    * len);
    }

    resize_carray(&tau->T_in, tau->start, len);
    resize_carray(&tau->T_out, tau->start, len);
    if (tau->use_fwd)
        resize_carray(&tau->T_fwd, tau->start, len);

    resize_array(&tau->rho_km_vals, tau->start, len);
    resize_array(&tau->F_km_vals, tau->start, len);
    resize_array(&tau->phi_rad_vals, tau->start, len);
    resize_array(&tau->k_vals, tau->start, len);
    resize_array(&tau->f_sky_hz_vals, tau->start, len);
    resize_array(&tau->rho_dot_kms_vals, tau->start, len);
    resize_array(&tau->raw_tau_threshold_vals, tau->start, len);
    resize_array(&tau->B_rad_vals, tau->start, len);
    resize_array(&tau->D_km_vals, tau->start, len);
    resize_array(&tau->w_km_vals, tau->start, len);
    resize_array(&tau->t_oet_spm_vals, tau->start, len);
    resize_array(&tau->t_ret_spm_vals, tau->start, len);
    resize_array(&tau->t_set_spm_vals, tau->start, len);
    resize_array(&tau->rho_corr_pole_km_vals, tau->start, len);
    resize_array(&tau->rho_corr_timing_km_vals, tau->start, len);
    resize_array(&tau->phi_rl_rad_vals, tau->start, len);
    resize_array(&tau->p_norm_vals, tau->start, len);
    resize_array(&tau->phase_rad_vals, tau->start, len);

    factor = log(tau->dx_km / tau->res);

    for (n = 0; n < len; ++n)
    {
        tau->power_vals[n] = tmpl_CDouble_Abs_Squared(tau->T_out[n]);
        tau->phase_vals[n] = tmpl_CDouble_Argument(tau->T_out[n]);

        mu = sin(fabs(tau->B_rad_vals[n]));
        tau->tau_vals[n] = -mu*log(tau->power_vals[n]);
        tau->tau_threshold_vals[n] = tau->raw_tau_threshold_vals[n] - factor*mu;

        if (tau->use_fwd)
        {
            tau->p_norm_fwd_vals[n] = tmpl_CDouble_Abs_Squared(tau->T_fwd[n]);

            /*  We multiplied by (1+i)/F instead of (1-i)/F, which is what    *
             *  the forward transformation wants. (1-i)/(1+i) = -i, which     *
             *  translates to a translation in phase of 3 pi / 2, or -pi/2.   */
            tau->phase_fwd_vals[n] = tmpl_CDouble_Argument(tau->T_fwd[n]) -
                                      tmpl_Pi_By_Two;
            tau->tau_fwd_vals[n] = -mu*log(tau->p_norm_fwd_vals[n]);
        }
    }

    tau->arr_size = len;
}

#endif
