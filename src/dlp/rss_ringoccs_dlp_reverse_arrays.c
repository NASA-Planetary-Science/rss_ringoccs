/*  tmpl_Double_Array_Reverse declared here, reverse the order of arrays.     */
#include <libtmpl/include/tmpl_math.h>

/*  Header file with the DLPObj typedef.                                      */
#include <rss_ringoccs/include/types/rss_ringoccs_dlpobj.h>

/*  size_t typedef given here.                                                */
#include <stddef.h>

/*  puts function found here, used for printing a status message if requested.*/
#include <stdio.h>

/*  Forward declaration / function prototype.                                 */
extern void rssringoccs_DLP_Reverse_Occultation(rssringoccs_DLPObj * const dlp);

/*  Reverses all of the arrays in a DLP object.                               */
void rssringoccs_DLP_Reverse_Occultation(rssringoccs_DLPObj * const dlp)
{
    /*  Variable for indexing the rho_dot_kms_vals arrays.                    */
    size_t n;

    /*  If the input is NULL there is nothing to be done.                     */
    if (!dlp)
        return;

    /*  Do not attempt to inspect the data if an error has already occurred.  */
    if (dlp->error_occurred)
        return;

    /*  Print a status message if the user requested one.                     */
    if (dlp->verbose)
        puts("\r\tDLP: Reversing the order of the DLP arrays...");

    /*  Reverse the arrays one-by-one.                                        */
    tmpl_Double_Array_Reverse(dlp->rho_km_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->phi_deg_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->B_deg_vals,dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->D_km_vals,dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->rho_dot_kms_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->t_oet_spm_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->t_ret_spm_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->t_set_spm_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->phi_rl_deg_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->rx_km_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->ry_km_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->rz_km_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->rho_corr_pole_km_vals, dlp->arr_size);
    tmpl_Double_Array_Reverse(dlp->rho_corr_timing_km_vals, dlp->arr_size);

    /*  We've swapped the order, flip the direction of rho_kms_vals.          */
    for(n = 0; n < dlp->arr_size; ++n)
        dlp->rho_dot_kms_vals[n] = -dlp->rho_dot_kms_vals[n];

    /*  Similarly, the displacement now moves in the opposite direction.      */
    dlp->dx_km = -dlp->dx_km;
}
/*  End of rssringoccs_DLP_Reverse_Occultation.                               */
