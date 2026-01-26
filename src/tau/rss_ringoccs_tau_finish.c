
#include <libtmpl/include/types/tmpl_complex_double.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <stdlib.h>

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
