/*  Booleans provided by this library.                                        */
#include <libtmpl/include/tmpl.h>

/*  Header file with the Tau definition and function prototype.               */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  NULL pointer is given here.                                               */
#include <stddef.h>

/*  Function for initializing all members in a tau object to NULL.            */
void rssringoccs_Tau_Init(rssringoccs_TAUObj *tau)
{
    /*  If the input is a NULL pointer there is nothing to be done. Abort.    */
    if (!tau)
        return;

    /*  Fresh tau object, so set the error values to their zero values.       */
    tau->error_occurred = tmpl_False;
    tau->error_message = NULL;

    /*  Initialize all pointers to NULL. This prevents things like double     *
     *  double free's, leaking memory by calling malloc twice, etc.           */
    tau->dlp = NULL;
    tau->T_in = NULL;
    tau->T_out = NULL;
    tau->T_fwd = NULL;
    tau->F_km_vals = NULL;
    tau->k_vals = NULL;
    tau->w_km_vals = NULL;
    tau->tau_threshold_vals = NULL;

    /*  Set the indexing variables to be zero as well.                        */
    tau->start = 0;
    tau->n_used = 0;

    /*  Set the remaining variables to their defaults.                        */
    rssringoccs_Tau_Set_Default_Values(tau);
}
/*  End of rssringoccs_Tau_Init.                                              */
