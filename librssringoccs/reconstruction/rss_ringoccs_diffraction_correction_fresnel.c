#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_bool.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <rss_ringoccs/include/rss_ringoccs_special_functions.h>
#include <rss_ringoccs/include/rss_ringoccs_error_check.h>

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
void DiffractionCorrectionFresnel(rssringoccs_TAUObj *dlp)
{
    /*  If everything executes smoothly, status should remain at zero.        */
    dlp->status = 0;

    /*  Check that the pointers to the data are not NULL.                     */
    if (check_dlp_data(dlp) == false)
    {
        /*  One of the variables has null data, return to calling function.   *
         *  Note, it is the callers responsibility to handle dlp->status      *
         *  accordingly and raise the appropriate errors.                     */
        dlp->status = 1;
        return;
    }

    /*  m and n used for indexing, nw_pts is number of points in window.      */
    long m, n, nw_pts, center;

    /*  w_init is window width (km), dx and two_dx are sample spacing (km).   */
    double w_init, dx, two_dx;

    /*  Pointers for the independent variable and the window function.        */
    double *x_arr;
    double *w_func;

    /*  If forward transform is selected, we'll need to negate x_arr. The     *
     *  following is a quick way to multiply by +/- 1. (dlp->use_fwd == 0) is *
     *  a Boolean which will be 1 if the statement is true and 0 if false.    *
     *  Similarly for (dlp->use_fwd == 1). Hence this returns +/- 1.          */
    int fwd_factor = (dlp->use_fwd == 0) - (dlp->use_fwd == 1);

    /*  Declare the window function pointer and allocate memory for it. The   *
     *  type rss_ringoccs_window_func was declared at the start of this file. *
     *  Be sure to free this at the end!                                      */
    rss_ringoccs_window_func fw;
    fw = malloc(sizeof(fw));

    /*  Cast the selected window type to the fw pointer.                      */
    select_window_func(fw, dlp);

    /*  select_window_func sets dlp->status to 5 if it could not successfully *
     *  parse dlp->wtype. As stated, its our responsibility to check this     *
     *  whenever using the function. Check and return to caller if it failed. */
    if (dlp->error_occured)
        return;

    /*  Cast FresT to the appropriate function.                               */
    if (dlp->use_norm)
        FresT = &Fresnel_Transform_Norm_Double;
    else
        FresT = &Fresnel_Transform_Double;

    /*  Retrieve the starting point from the DLPObj instance.                 */
    center = dlp->start;

    /*  Compute some extra necessary variables. At this point it is assumed   *
     *  that dlp->w_km_vals, dlp->rho_km_vals, and others are pointers, most  *
     *  likely created with malloc or calloc, that point to a memory block    *
     *  that is dlp->arr_size in size, that dlp->start >= 0, and that         *
     *  dlp->start+dlp->n_used <= dlp->arr_size. No error checks for this     *
     *  are performed here, but rather the caller of this function has that   *
     *  responsibility. Such checks are performed in the                      *
     *  DiffractionCorrection Python class, so if you're only using that then *
     *  there's no problem. If not, this next step may cause a segmentation   *
     *  fault.                                                                */
    w_init = dlp->w_km_vals[center];

    /*  It is also assumed these pointers have at least two elements of       *
     *  doubles being pointed to. Again, DiffractionCorrection checks this.   *
     *  Hence dlp->rho_km_vals[center] and dlp->rho_km_vals[center+1] should  *
     *  be initialized. If not, you will get a segmentation fault.            */
    dx = dlp->rho_km_vals[center+1] - dlp->rho_km_vals[center];

    /*  And now, some more variables we'll be using.                          */
    two_dx  = 2.0*dx;
    nw_pts  = ((long)(w_init / two_dx))+1;

    /* Check to ensure you have enough data to the left.                      */
    if (!check_data_range(dlp, two_dx))
    {
        /*  One of the points has too large of a window width to process.     *
         *  Returning with error message. Again, it is the callers job to     *
         *  check dlp->status for errors.                                     */
        dlp->status = 2;
        return;
    }

    /*  Reserve some memory for two arrays, the ring radius and the window    *
     *  function. This will need to be reallocated later if the window width  *
     *  changes by more than two_dx. Per C99 recommendations, we do not cast  *
     *  malloc since void will safely be promoted to a double pointer.        */
    x_arr  = malloc(sizeof(*x_arr)  * nw_pts);
    w_func = malloc(sizeof(*w_func) * nw_pts);

    /*  Check that malloc was successfull then pass the x_arr array           *
     *  (ring radius) to the void function reset_window. This alters x_arr so *
     *  that it's values range from -W/2 to zero, W begin the window width.   */
    if (!(x_arr) || !(w_func))
    {
        /*  Malloc failed, return to calling function.                        */
        dlp->status = 3;
        return;
    }
    else
        reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

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
    for (m=0; m<=dlp->n_used; ++m)
    {
        /*  If the window width has deviated more the 2*dx, reset variables.  *
         *  fabs is the absolute value function for double precision          *
         *  variables and is defined in the built-in math.h.                  */
        if (fabs(w_init - dlp->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init = dlp->w_km_vals[center];
            nw_pts = ((long)(w_init / two_dx))+1;

            /*  Reallocate memory, since the sizes of the arrays changed.     */
            w_func = (double *)realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = (double *)realloc(x_arr, sizeof(double)*nw_pts);

            /*  Reset the x_arr array to range between -W/2 and zero.         */
            reset_window(x_arr, w_func, dx, w_init, nw_pts, fw);

            /* Compute Window Functions, and compute pi/2 * x^2               */
            for(n=0; n<nw_pts; ++n)
            {
                x_arr[n] *= rssringoccs_Pi_By_Two*x_arr[j];

                /*  Again, if forward calculation is set, negate x_arr.       */
                x_arr[n] *= fwd_factor;
            }
        }

        /*  Compute the Fresnel Transform about the current point.            */
        dlp->T_out[m] = FresT(x_arr, dlp->T_in, w_func,
                              dlp->F_km_vals[center], dx, nw_pts, center);

        /*  Move the pointers to the next point.                              */
        center += 1;
    }

    /*  Free the variables allocated by malloc.                               */
    free(x_arr);
    free(w_func);
    free(fw);
}
