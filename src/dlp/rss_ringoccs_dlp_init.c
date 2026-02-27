/*  NULL pointer is given here.                                               */
#include <stddef.h>

/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  Function for initializing all members in a tau object to NULL.            */
void rssringoccs_Tau_Init(rssringoccs_TAUObj *tau)
{
    /*  Variable for zero cast to type "size_t".                              */
    const size_t zero = (size_t)0;

    /*  If the input is a NULL pointer there is nothing to be done. Abort.    */
    if (!tau)
        return;

    /*  Fresh tau object, so set the error values to their zero values.       */
    tau->error_occurred = tmpl_False;
    tau->error_message = NULL;

    /*  Initialize all pointers to NULL. This prevents things like double     *
     *  double free's, leaking memory by calling malloc twice, etc.           */
    tau->T_in = NULL;
    tau->T_out = NULL;
    tau->T_fwd = NULL;
    tau->rho_km_vals = NULL;
    tau->F_km_vals = NULL;
    tau->phi_deg_vals = NULL;
    tau->k_vals = NULL;
    tau->rho_dot_kms_vals = NULL;
    tau->B_deg_vals = NULL;
    tau->D_km_vals = NULL;
    tau->w_km_vals = NULL;
    tau->t_oet_spm_vals = NULL;
    tau->t_ret_spm_vals = NULL;
    tau->t_set_spm_vals = NULL;
    tau->rho_corr_pole_km_vals = NULL;
    tau->rho_corr_timing_km_vals = NULL;
    tau->tau_threshold_vals = NULL;
    tau->phi_rl_deg_vals = NULL;
    tau->rx_km_vals = NULL;
    tau->ry_km_vals = NULL;
    tau->rz_km_vals = NULL;

    /*  Set the indexing variables to be zero as well.                        */
    tau->arr_size = zero;
    tau->start = zero;
    tau->n_used = zero;

    /*  Set the remaining variables to their defaults.                        */
    rssringoccs_Tau_Set_Default_Values(tau);
}
/*  End of rssringoccs_Tau_Init.                                              */
