#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <stdlib.h>

void rssringoccs_Tau_Get_Window_Width(rssringoccs_TAUObj* tau)
{
    /*  Declare long pointer-to-pointer which stores the indices where        *
     *  F_km_vals is non-zero in the first slot (Prange[0]), and the size of  *
     *  this array in the second (*Prange[1]).                                */
    unsigned long **Prange, **wrange;
    unsigned long *Prange_Index, *wrange_Index;
    unsigned long Prange_Size, wrange_Size;
    unsigned long n;
    double w_fac, omega, min_val, max_val;
    double *alpha, *P_vals, *rho_legal;

    if (tau == NULL)
        return;
    else if (tau->error_occurred)
        return;

    if (tau->w_km_vals != NULL)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rtau->w_km_vals is not NULL. You have either already processed\n"
            "\rthe window width or did not create the rssringoccs_TAUObj\n"
            "\rpointer with rssringoccs_Create_TAUObj. Returning with error.\n"
        );
        return;
    }

    /*  Use calloc to both allocate memory for tau.w_km_vals (like malloc)    *
     *  and initialize the data to zero (unlike malloc). This is similar to   *
     *  numpy.zeros(tau.arr_size) in Python.                                  */
    tau->w_km_vals = calloc(tau->arr_size, sizeof(double));

    if (tau->bfac)
    {
        w_fac = tau->normeq;
        alpha  = malloc(sizeof(*alpha) * tau->arr_size);
        P_vals = malloc(sizeof(*P_vals) * tau->arr_size);

        for(n=0; n<tau->arr_size; ++n)
        {
            omega     = rssringoccs_Two_Pi * tau->f_sky_hz_vals[n];
            alpha[n]  = omega * tau->sigma;
            alpha[n] *= alpha[n] * 0.5 / tau->rho_dot_kms_vals[n];
            P_vals[n] = tau->res/(alpha[n]*tau->F_km_vals[n]*tau->F_km_vals[n]);
        }

        Prange = rssringoccs_Where_Greater_Double(P_vals, tau->arr_size, 1.0);
        Prange_Index = Prange[0];
        Prange_Size  = *Prange[1];

        for (n=0; n<Prange_Size; ++n)
            tau->w_km_vals[Prange_Index[n]] = w_fac *
                rssringoccs_Double_Resolution_Inverse(P_vals[Prange_Index[n]]) /
                alpha[n];

        free(P_vals);
        free(alpha);
    }
    else
    {
        w_fac = tau->normeq/tau->res;

        for (n=0; n < tau->arr_size; ++n)
            tau->w_km_vals[n] = 2.0*tau->F_km_vals[n]*tau->F_km_vals[n]*w_fac;

        Prange = rssringoccs_Where_Greater_Double(tau->F_km_vals,
                                                  tau->arr_size, 0.0);
        Prange_Index = Prange[0];
        Prange_Size  = *Prange[1];
    }

    rho_legal = malloc(sizeof(*rho_legal) * Prange_Size);

    for(n=0; n<Prange_Size; ++n)
        rho_legal[n] = tau->rho_km_vals[Prange_Index[n]] +
                       0.5*tau->w_km_vals[Prange_Index[n]];

    min_val = rssringoccs_Min_Double(rho_legal, Prange_Size);
    max_val = rssringoccs_Max_Double(rho_legal, Prange_Size);
    wrange = rssringoccs_Where_LesserGreater_Double(tau->rho_km_vals,
                                                    tau->arr_size,
                                                    min_val, max_val);
    wrange_Index = wrange[0];
    wrange_Size = *wrange[1];
    free(Prange_Index);
    free(Prange);
    free(rho_legal);

    if (wrange_Size == 0)
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\r\tThe window width is too large to reconstruct anything.\n"
        );
        return;
    }
    else if (max_val < tau->rng_req[0])
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\r\tMinimum requested range is greater than available data.\n"
        );
        free(wrange_Index);
        free(wrange);
        return;
    }
    else if (min_val > tau->rng_req[1])
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\r\tMinimum requested range is greater than available data.\n"
        );
        free(wrange_Index);
        free(wrange);
        return;
    }
    else
    {
        free(wrange_Index);
        free(wrange);
        tau->start = wrange_Index[0];
        tau->n_used = wrange_Index[wrange_Size-1] - tau->start;
    }

    if (tau->start > tau->arr_size){
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rStarting index (start) is greater than the size of the array.\n"
        );
        return;
    }
    else if (tau->start + tau->n_used > tau->arr_size){
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rFinal index (start+n_used) is greater than size of array.\n"
        );
        return;
    }
}
