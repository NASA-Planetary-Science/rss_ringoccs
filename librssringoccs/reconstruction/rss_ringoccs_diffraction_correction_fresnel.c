#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionFresnel                                          *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the classic Fresnel quadratic     *
 *      approximation to the Fresnel kernel.                                  *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          rss_ringoccs_diffraction_correction.h. This contains all of the   *
 *          necessary data for diffraction correction, including the geometry *
 *          of the occultation and actual power and phase data.               *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction corrected data.                   *
 *  Notes:                                                                    *
 *      1.) This code uses the Fresnel approximation which has been known to  *
 *          fail for several different occultations, especially ones of very  *
 *          low angle (small B values). Take this into consideration when     *
 *          performing any analysis.                                          *
 *      2.) While this may be inaccurate for certain occultations, it is      *
 *          immensely fast, capable of processing the entire Rev007 E         *
 *          occultation accurately in less than a second at 1km resolution.   *
 ******************************************************************************/
void rssringoccs_Diffraction_Correction_Fresnel(rssringoccs_TAUObj *tau)
{
    /*  m and n used for indexing, nw_pts is number of points in window.      */
    unsigned long m, n, nw_pts, center;

    /*  w_init is window width (km), dx and two_dx are sample spacing (km).   */
    double w_init, dx, two_dx;

    /*  Pointers for the independent variable and the window function.        */
    double *x_arr;
    double *w_func;
    double fwd_factor;

    /*  Declare the window function pointer and allocate memory for it. The   *
     *  type rssringoccs_window_func was declared at the start of this file.  *
     *  Be sure to free this at the end!                                      */
    rssringoccs_window_func fw = tau->window_func;

    /*  This should remain at false.                                          */
    tau->error_occurred = rssringoccs_False;

    if (tau->use_fwd)
        fwd_factor = -1.0;
    else
        fwd_factor = 1.0;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Check_Tau_Data(tau);
    if (tau->error_occurred)
        return;

    /*  Retrieve the starting point from the TAUObj instance.                 */
    center = tau->start;

    /*  Compute some extra necessary variables. At this point it is assumed   *
     *  that tau->w_km_vals, tau->rho_km_vals, and others are pointers, most  *
     *  likely created with malloc or calloc, that point to a memory block    *
     *  that is tau->arr_size in size, that tau->start >= 0, and that         *
     *  tau->start+tau->n_used <= tau->arr_size. No error checks for this     *
     *  are performed here, but rather the caller of this function has that   *
     *  responsibility. Such checks are performed in the                      *
     *  DiffractionCorrection Python class, so if you're only using that then *
     *  there's no problem. If not, this next step may cause a segmentation   *
     *  fault.                                                                */
    w_init = tau->w_km_vals[center];

    /*  It is also assumed these pointers have at least two elements of       *
     *  doubles being pointed to. Again, DiffractionCorrection checks this.   *
     *  Hence dlp->rho_km_vals[center] and dlp->rho_km_vals[center+1] should  *
     *  be initialized. If not, you will get a segmentation fault.            */
    dx = tau->rho_km_vals[center+1] - tau->rho_km_vals[center];

    /*  And now, some more variables we'll be using.                          */
    two_dx = 2.0*dx;
    nw_pts = ((long)(w_init / two_dx)) + 1;

    /* Check to ensure you have enough data to the left.                      */
    rssringoccs_Check_Tau_Data_Range(tau);
    if (tau->error_occurred)
        return;

    /*  Reserve some memory for two arrays, the ring radius and the window    *
     *  function. This will need to be reallocated later if the window width  *
     *  changes by more than two_dx. Per C99 recommendations, we do not cast  *
     *  malloc since void will safely be promoted to a double pointer.        */
    x_arr  = malloc(sizeof(*x_arr)  * nw_pts);
    w_func = malloc(sizeof(*w_func) * nw_pts);

    /*  Check that malloc was successfull then pass the x_arr array           *
     *  (ring radius) to the void function reset_window. This alters x_arr so *
     *  that it's values range from -W/2 to zero, W begin the window width.   */
    if (!(x_arr))
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tDiffractionCorrectionFresnel\n\n"
            "\rMalloc failed and returned NULL for x_arr. Returning.\n\n"
        );
        return;
    }

    if (!(w_func))
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tDiffractionCorrectionFresnel\n\n"
            "\rMalloc failed and returned NULL for w_func. Returning.\n\n"
        );
        return;
    }

    rssringoccs_Tau_Reset_Window(x_arr, w_func, dx, w_init, nw_pts, fw);

    /* Compute Window Functions, and compute pi/2 * x^2                       */
    for(m=0; m<nw_pts; ++m)
    {
        /*  The independent variable is pi/2 * ((rho-rho0)/F)^2. Compute      *
         *  part of this. The 1/F^2 part is introduced later.                 */
        x_arr[m] *= rssringoccs_Pi_By_Two*x_arr[m];

        /*  Use the fwd_factor to computer forward or inverse transform.      */
        x_arr[m] *= fwd_factor;
    }

    /*  Compute the Fresnel transform across the input data.                  */
    if (tau->use_norm)
    {
        /*  Compute the Fresnel transform across the input data.              */
        for (m=0; m<=tau->n_used; ++m)
        {
            /*  If the window width has deviated more the 2*dx, reset         *
             *  variables. fabs is the absolute value function for double     *
             *  precision variables and is defined in the built-in math.h.    */
            if (fabs(w_init - tau->w_km_vals[center]) >= two_dx)
            {
                /* Reset w_init and recompute window function.                */
                w_init = tau->w_km_vals[center];
                nw_pts = ((long)(w_init / two_dx))+1;

                /*  Reallocate memory, since the sizes of the arrays changed. */
                w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
                x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);

                /*  Reset the x_arr array to range between -W/2 and zero.     */
                rssringoccs_Tau_Reset_Window(x_arr, w_func, dx,
                                             w_init, nw_pts, fw);

                /* Compute Window Functions, and compute pi/2 * x^2.          */
                for(n=0; n<nw_pts; ++n)
                {
                    x_arr[n] *= rssringoccs_Pi_By_Two*x_arr[n];

                    /*  Again, if forward calculation is set, negate x_arr.   */
                    x_arr[n] *= fwd_factor;
                }
            }

            /*  Compute the Fresnel Transform about the current point.        */
            tau->T_out[m] = Fresnel_Transform_Norm_Double(
                x_arr, tau->T_in, w_func,
                tau->F_km_vals[center], nw_pts, center
            );

            /*  Move the pointers to the next point.                          */
            center += 1;
        }
    }
    else
    {
        for (m=0; m<=tau->n_used; ++m)
        {
            /*  If the window width has deviated more the 2*dx, reset         *
             *  variables. fabs is the absolute value function for double     *
             *  precision variables and is defined in the built-in math.h.    */
            if (fabs(w_init - tau->w_km_vals[center]) >= two_dx)
            {
                /* Reset w_init and recompute window function.                */
                w_init = tau->w_km_vals[center];
                nw_pts = ((long)(w_init / two_dx))+1;

                /*  Reallocate memory, since the sizes of the arrays changed. */
                w_func = realloc(w_func, sizeof(*w_func) * nw_pts);
                x_arr  = realloc(x_arr,  sizeof(*x_arr)  * nw_pts);

                /*  Reset the x_arr array to range between -W/2 and zero.     */
                rssringoccs_Tau_Reset_Window(x_arr, w_func, dx,
                                             w_init, nw_pts, fw);

                /* Compute Window Functions, and compute pi/2 * x^2.          */
                for(n=0; n<nw_pts; ++n)
                {
                    x_arr[n] *= rssringoccs_Pi_By_Two*x_arr[n];

                    /*  Again, if forward calculation is set, negate x_arr.   */
                    x_arr[n] *= fwd_factor;
                }
            }

            /*  Compute the Fresnel Transform about the current point.        */
            tau->T_out[m] = Fresnel_Transform_Double(x_arr, tau->T_in, w_func,
                                                    tau->F_km_vals[center], dx,
                                                    nw_pts, center);

            /*  Move the pointers to the next point.                          */
            center += 1;
        }
    }

    /*  Free the variables allocated by malloc.                               */
    free(x_arr);
    free(w_func);
}
