#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_optics.h>
#include <libtmpl/include/tmpl_string.h>
#include <libtmpl/include/tmpl_special_functions_real.h>
#include <libtmpl/include/tmpl_where.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stdlib.h>
#include <stdio.h>

void rssringoccs_Tau_Get_Window_Width(rssringoccs_TAUObj* tau)
{
    /*  Declare long pointer-to-pointer which stores the indices where        *
     *  F_km_vals is non-zero in the first slot (Prange[0]), and the size of  *
     *  this array in the second (*Prange[1]).                                */
    size_t **Prange, **wrange;
    size_t *Prange_Index, *wrange_Index;
    size_t Prange_Size, wrange_Size;
    size_t n;
    double w_fac, omega, F;
    double *alpha, *P_vals, *rho_legal;

    if (tau == NULL)
        return;

    if (tau->error_occurred)
        return;

    if (tau->rng_list[0] <= tau->rho_km_vals[0])
    {
        tau->rng_list[0] = tau->rho_km_vals[0];
        tau->start = 0;
    }
    else if (tau->rng_list[0] > tau->rho_km_vals[tau->arr_size-1])
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = malloc(sizeof(*tau->error_message)*512);
        if (tau->error_message == NULL)
            return;

        sprintf(
            tau->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rMinimum requested radius is greater than the largest\n"
            "\rradius available in the data.\n"
            "\r\tRequested Radii: %f %f\n"
            "\r\tAvailable Radii: %f %f\n",
            tau->rng_list[0], tau->rng_list[1],
            tau->rho_km_vals[0], tau->rho_km_vals[tau->arr_size-1]
        );
        return;
    }
    else
    {
        n = 0;
        while (tau->rho_km_vals[n] < tau->rng_list[0])
            n++;

        tau->rng_list[0] = tau->rho_km_vals[n];
        tau->start = n;
    }

    if (tau->rng_list[1] >= tau->rho_km_vals[tau->arr_size-1])
    {
        tau->rng_list[1] = tau->rho_km_vals[tau->arr_size-1];
        tau->n_used = tau->arr_size - tau->start;
    }
    else if (tau->rng_list[1] < tau->rho_km_vals[0])
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = malloc(sizeof(*tau->error_message)*512);
        if (tau->error_message == NULL)
            return;

        sprintf(
            tau->error_message,
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rMinimum requested radius is greater than the largest\n"
            "\rradius available in the data.\n"
            "\r\tRequested Radii: %f %f\n"
            "\r\tAvailable Radii: %f %f\n",
            tau->rng_list[0], tau->rng_list[1],
            tau->rho_km_vals[0], tau->rho_km_vals[tau->arr_size-1]
        );
        return;
    }
    else
    {
        n = tau->start;
        while (tau->rho_km_vals[n] <= tau->rng_list[1])
            n++;

        tau->rng_list[1] = tau->rho_km_vals[n];
        tau->n_used = n - tau->start;
    }

    if (tau->w_km_vals != NULL)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
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
        alpha  = malloc(sizeof(*alpha) * tau->n_used);
        P_vals = malloc(sizeof(*P_vals) * tau->n_used);

        for(n = 0; n < tau->n_used; ++n)
        {
            F = tau->F_km_vals[n + tau->start];
            omega = tmpl_Speed_Of_Light_KMS * tau->k_vals[n+tau->start];
            alpha[n] = omega * tau->sigma;
            alpha[n] *= alpha[n] * 0.5 / tau->rho_dot_kms_vals[n];
            P_vals[n] = tau->res/(alpha[n]*F*F);
        }

        Prange = tmpl_Where_Greater_Double(P_vals, tau->n_used, 1.0);
        Prange_Index = Prange[0];
        Prange_Size  = *Prange[1];

        for (n=0; n<Prange_Size; ++n)
            tau->w_km_vals[tau->start + Prange_Index[n]] = w_fac *
                tmpl_Double_Resolution_Inverse(P_vals[Prange_Index[n]]) /
                alpha[n];

        free(P_vals);
        free(alpha);
        free(Prange_Index);
        free(Prange);
    }
    else
    {
        w_fac = tau->normeq/tau->res;

        for (n=0; n < tau->n_used; ++n)
        {
            F = tau->F_km_vals[n + tau->start];
            tau->w_km_vals[n+tau->start] = 2.0*F*F*w_fac;
        }
    }

    rho_legal = malloc(sizeof(*rho_legal) * tau->n_used);

    for(n=0; n<tau->n_used; ++n)
        rho_legal[n] = tau->rho_km_vals[tau->start + n] -
                       0.5*tau->w_km_vals[tau->start + n];

    wrange = tmpl_Where_Greater_Double(rho_legal, tau->n_used,
                                       tau->rho_km_vals[0]);
    wrange_Index = wrange[0];
    wrange_Size = *wrange[1];

    if (wrange_Size == 0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rThe window width is too large to reconstruct anything.\n"
            "\rThere is not enough data to the left of any point. That is,\n"
            "\rrho_km_vals[n] - w_km_vals[n]/2 is less than rho_km_vals[0]\n"
            "\rfor all n. Returning with error.\n"
        );
        free(wrange_Index);
        free(wrange);
        return;
    }
    else if (wrange_Index[0] > 0)
    {
        tau->start  += wrange_Index[0];
        tau->n_used -= wrange_Index[0];
    }

    free(wrange_Index);
    free(wrange);

    for(n=0; n<tau->n_used; ++n)
        rho_legal[n] = tau->rho_km_vals[tau->start + n] +
                       0.5*tau->w_km_vals[tau->start + n];

    wrange = tmpl_Where_Lesser_Double(rho_legal, tau->n_used,
                                      tau->rho_km_vals[tau->arr_size-1]);

    wrange_Index = wrange[0];
    wrange_Size = *wrange[1];
    free(rho_legal);

    if (wrange_Size == 0)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rThe window width is too large to reconstruct anything.\n"
            "\rThere is not enough data to the left of any point. That is,\n"
            "\rrho_km_vals[n] + w_km_vals[n]/2 is less than rho_km_vals[0]\n"
            "\rfor all n. Returning with error.\n"
        );
        return;
    }
    else if (wrange_Index[wrange_Size-1] < tau->n_used)
        tau->n_used = wrange_Size-1;

    free(wrange_Index);
    free(wrange);

    if (tau->start > tau->arr_size)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = malloc(sizeof(*tau->error_message)*512);
        if (tau->error_message == NULL)
            return;

        sprintf(
            tau->error_message,
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rStarting index (start) is greater than the size of the array.\n"
            "\r\ttau->start:    %lu\n"
            "\r\ttau->n_used:   %lu\n"
            "\r\tstart+n_used:  %lu\n"
            "\r\ttau->arr_size: %lu\n",
            tau->start, tau->n_used, tau->start+tau->n_used, tau->arr_size
        );
        return;
    }
    else if ((tau->start + tau->n_used) > tau->arr_size)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = malloc(sizeof(*tau->error_message)*512);
        if (tau->error_message == NULL)
            return;

        sprintf(
            tau->error_message,
            "\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Get_Window_Width\n\n"
            "\rFinal index (start+n_used) is greater than size of array.\n"
            "\r\ttau->start:    %lu\n"
            "\r\ttau->n_used:   %lu\n"
            "\r\tstart+n_used:  %lu\n"
            "\r\ttau->arr_size: %lu\n",
            tau->start, tau->n_used, tau->start+tau->n_used, tau->arr_size
        );
        return;
    }
}
