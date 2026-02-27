/*  NULL pointer is given here.                                               */
#include <stddef.h>

/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the DLP definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_dlp.h>

/*  Function for initializing all members in a dlp object to NULL.            */
void rssringoccs_DLP_Init(rssringoccs_DLPObj * const dlp)
{
    /*  If the input is a NULL pointer there is nothing to be done. Abort.    */
    if (!dlp)
        return;

    /*  Initialize all pointers to NULL. This prevents things like double     *
     *  double free's, leaking memory by calling malloc twice, etc.           */
    dlp->rho_km_vals = NULL;
    dlp->phi_deg_vals = NULL;
    dlp->B_deg_vals = NULL;
    dlp->D_km_vals = NULL;
    dlp->f_sky_hz_vals = NULL;
    dlp->rho_dot_kms_vals = NULL;
    dlp->t_oet_spm_vals = NULL;
    dlp->t_ret_spm_vals = NULL;
    dlp->t_set_spm_vals = NULL;
    dlp->rho_corr_pole_km_vals = NULL;
    dlp->rho_corr_timing_km_vals = NULL;
    dlp->phi_rl_deg_vals = NULL;
    dlp->p_norm_vals = NULL;
    dlp->phase_deg_vals = NULL;
    dlp->raw_tau_threshold_vals = NULL;
    dlp->rx_km_vals = NULL;
    dlp->ry_km_vals = NULL;
    dlp->rz_km_vals = NULL;

    dlp->dx_km = 0.0;

    /*  Set the indexing variables to be zero as well.                        */
    dlp->arr_size = 0;
    dlp->reference_count = 0;

    /*  Fresh DLP object, set the error variables to their zero values.       */
    dlp->error_occurred = tmpl_False;
    dlp->error_message = NULL;
}
/*  End of rssringoccs_DLP_Init.                                              */
