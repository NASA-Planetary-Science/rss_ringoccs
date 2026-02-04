#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <libtmpl/include/tmpl_math.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Performs the Newton transform on a Tau object.                            */
void rssringoccs_Diffraction_Correction_Newton(rssringoccs_TAUObj * const tau)
{
    /*  Variable for indexing the inner most loop.                            */
    size_t n;

    /*  Check that tau isn't NULL before trying to access its members.        */
    if (!tau)
        return;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Core_Data(tau);

    /* Check to ensure you have enough data to process the transform.         */
    rssringoccs_Tau_Check_Data_Range(tau);

    /*  The previous functions set the error_occurred Boolean to true on      *
     *  error. Do not proceed if an error occurred, check for this.           */
    if (tau->error_occurred)
        return;

    {
        /*  The rho_km_vals array is required to have equally spaced samples. *
         *  The displacement between bins can be computed from the difference *
         *  of any two successive points.                                     */
        const double dx = tau->rho_km_vals[1] - tau->rho_km_vals[0];
        const double two_dx = 2.0 * dx;

        /*  When the window width is large, or the sample spacing is very fine,   *
         *  repeated calls to realloc, which resizes the window array, start to   *
         *  bottleneck the computation. It is faster to compute the maximum width *
         *  required for the entire computation and then make one call to malloc. *
         *  These variables are for the largest window width required.            */
        const double * const w_array = tau->w_km_vals + tau->start;
        const double w_max = tmpl_Double_Array_Max(w_array, tau->n_used);
        const size_t nw_max = (TMPL_CAST(w_max / two_dx, size_t) << 1) + 1;

        /*  Set the correct function pointer.                                 */
        const rssringoccs_FresnelNewtonTransform newton_transform =
            rssringoccs_Tau_Select_Newton_Transform(tau);

        /*  The selector function checks if the psinum in the Tau object does *
         *  not represent a Newton transform. It will set the error_occurred  *
         *  Boolean should this take place. Check.                            */
        if (tau->error_occurred)
            return;

        /*  At the time of this writing, not every Newton transform has been  *
         *  implemented (there are many!). If the selector function returned  *
         *  NULL, the user requested a psitype that is not available.         */
        if (!newton_transform)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Diffraction_Correction_Newton\n\n"
                "\rrssringoccs_Tau_Select_Newton_Transform returned NULL.\n"
                "\rRequested psitype not implemented yet.\n";

            return;
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            /*  w_init is the initial window width. We do not need to repeatedly      *
             *  update the window function, we only need to recompute once the width  *
             *  grows by at least two bins (one on the left and one on the right).    *
             *  That is, once abs(w_km_vals[center] - w_init) > 2 * dx is true.       */
            double w_init = tau->w_km_vals[tau->start];

            /*  The windows must contain an odd number of samples, half to    *
             *  the left of the center, half to the right, and the center.    */
            size_t nw_pts = (TMPL_CAST(w_init / two_dx, size_t) << 1) + 1;


            /*  Pointer to a double array for the tapering function. The value        *
             *  w_func[n] corresponds to the evaluation of w(x[n]) at the nth point   *
             *  in the window. w is often some combination of Bessel functions, it is *
             *  beneficial to precompute the tapering function once and use that.     */
            double *w_func = TMPL_MALLOC(double, nw_max);

            /*  Malloc returns NULL on failure. Check for this.               */
            if (!w_func)
            {
                tau->error_occurred = tmpl_True;
                tau->error_message =
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\trssringoccs_Diffraction_Correction_Newton\n\n"
                    "\rmalloc returned NULL for w_func.\n";

                goto CLEANUP;
            }

            /*  Compute the values in the window function. Note, w_func almost surely *
             *  has more than nw_pts allocated to it, it has nw_max allocated. This   *
             *  function sets the first nw_pts in the array, the remaining points are *
             *  uninitialized. This makes no difference, all of the routines used     *
             *  take in the nw_pts variable as a parameter and do not read beyond     *
             *  that length in the w_func array.                                      */
            rssringoccs_Tau_Compute_Window(tau, w_func, nw_pts, tau->start);

            /*  Loop through all of the points in the processing range and run the    *
             *  transform (either forward or inverse) on the data.                    */
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (n = 0; n < tau->n_used; ++n)
            {
                /*  The window width varies with the Fresnel scale, which depends on  *
                 *  the geometry and hence varies across the data set. Compute how    *
                 *  much the width has changed so far.                                */
                const size_t center = tau->start + n;
                const double w_err = tmpl_Double_Abs(w_init - tau->w_km_vals[center]);

                /*  If the window width changes significantly, recompute w_func.      */
                if (w_err >= two_dx)
                {
                    /* Reset w_init and recompute window function.                    */
                    w_init = tau->w_km_vals[center];
                    nw_pts = (TMPL_CAST(w_init / two_dx, size_t) << 1) + 1;

                    /*  Compute the window function for the new window width.         */
                    rssringoccs_Tau_Compute_Window(tau, w_func, nw_pts, center);
                }

                /*  Compute the Newton transform about the current point.             */
                newton_transform(tau, w_func, nw_pts, center);
            }

            /*  Free variables allocated by malloc.                                   */
            CLEANUP:
                TMPL_FREE(w_func);
        }
    }
}
