#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  Performs the Newton transform using interpolating polynomials.            */
void
rssringoccs_Diffraction_Correction_Polynomial_Newton(
    rssringoccs_TAUObj * const tau
)
{
    /*  Variable for indexing the inner most loop.                            */
    size_t n;

    /*  Variable for the number of points in the window, and the index of the *
     *  central point in the window.                                          */
    size_t nw_pts, center;

    /*  w_init is the initial window width. We do not need to repeatedly      *
     *  update the window function, we only need to recompute once the width  *
     *  grows by at least two bins (one on the left and one on the right).    *
     *  That is, once abs(w_km_vals[center] - w_init) > 2 * dx is true.       */
    double w_init;

    /*  The distance, dx, between samples. 2 * dx is used frequently enough,  *
     *  we'll store this value in a variable as well.                         */
    double dx, two_dx;

    /*  Pointer to a double array for the tapering function. The value        *
     *  w_func[n] corresponds to the evaluation of w(x[n]) at the nth point   *
     *  in the window. w is often some combination of Bessel functions, it is *
     *  beneficial to precompute the tapering function once and use that.     */
    double *w_func = NULL;

    /*  Array for the dummy variable of integration, rho - rho0.              */
    double *x_arr = NULL;

    /*  Variable for the actual transform used. We'll select it later.        */
    rssringoccs_FresnelTransform newton_transform;

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

    /*  We have enough data to the left and right of the starting point to    *
     *  perform the transform. Initialize center to the starting index.       */
    center = tau->start;

    /*  We will not recompute the window function at each point. We update    *
     *  once this value differs from w_km_vals[n] by more than 2 * dx.        */
    w_init = tau->w_km_vals[center];

    /*  The rho_km_vals array is required to have equally spaced samples. The *
     *  displacement between bins can be computed from the difference of any  *
     *  two successive points.                                                */
    dx = tau->rho_km_vals[center + 1] - tau->rho_km_vals[center];
    two_dx = 2.0 * dx;

    /*  The windows must contain an odd number of samples, half to the left   *
     *  of the center, half to the right, and the center itself. From the     *
     *  symmetry of the interpolating polynomials, we do not need to          *
     *  explicitly store the right half or the center of the window in        *
     *  memory. nw_pts is hence the number of points in the left half only.   */
    nw_pts = TMPL_CAST(w_init / two_dx, size_t);

    /*  Set the correct function pointer.                                     */
    newton_transform = rssringoccs_Tau_Select_Newton_Interp_Transform(tau);

    /*  The selector function checks if the psinum in the Tau object does not *
     *  represent a Newton transform. It will set the error_occurred Boolean  *
     *  should this take place. Check.                                        */
    if (tau->error_occurred)
        return;

    /*  At the time of this writing, not every Newton transform has been      *
     *  implemented (there are many!). If the selector function returned NULL,*
     *  the user requested a psitype that is not currently available.         *
     *                                                                        *
     *      TODO:                                                             *
     *          Implement all psitypes and remove this error check.           */
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

    /*  Allocate memory for the independent variable and window function.     */
    x_arr = TMPL_MALLOC(double, nw_pts);
    w_func = TMPL_MALLOC(double, nw_pts);

    /*  Check to make sure malloc did not fail. It returns NULL on error.     */
    if (!w_func || ! x_arr)
    {
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Diffraction_Correction_Newton\n\n"
            "\rmalloc returned NULL for w_func.\n";

        /*  It is possible that malloc succeeded for one variable and not the *
         *  other. Since we initialized these pointers to NULL, if a pointer  *
         *  is not NULL, then malloc was successful in this case and we need  *
         *  to free it. The TMPL_FREE macro only calls free if the pointer is *
         *  not NULL. Use this on each variable.                              */
        TMPL_FREE(x_arr);
        TMPL_FREE(w_func);
        return;
    }

    /*  Initialize the window array and the independent variable.             */
    rssringoccs_Tau_Reset_Window(
        tau,                /*  Tau object containing the window function.    */
        x_arr,              /*  The independent variable, r[n]-r[center].     */
        w_func,             /*  The window array as a function of x_arr.      */
        nw_pts,             /*  Number of points in the x_arr array.          */
        center              /*  Index for the center of the window.           */
    );

    /*  Loop through all of the points in the processing range and run the    *
     *  transform (either forward or inverse) on the data.                    */
    for (n = 0; n < tau->n_used; ++n)
    {
        /*  Check if we need to resize the window. This happens once          *
         *  |w - w0| > 2 * dx occurs, where w is the current window width,    *
         *  and w0 is the value of w_init.                                    */
        const int resize = rssringoccs_Tau_Resize_Half_Window(
            tau, &x_arr, &w_func, w_init, two_dx, center
        );

        /*  If we did need a resize, there are a few things that could have   *
         *  gone wrong with memory reallocation, and a few things we'll need  *
         *  to reset if all of the memory management succeeded.               */
        if (resize)
        {
            /*  It is possible realloc failed during the call to              *
             *  rssringoccs_Tau_Resize_Half_Window. This function sets the    *
             *  error_occurred Boolean to True should this occur, and also    *
             *  frees the memory and sets the pointers to NULL. Simply abort  *
             *  the computation if any error occurred.                        */
            if (tau->error_occurred)
                return;

            /*  Reset the threshold value for the window width to the current *
             *  window size. We will update again once the required window    *
             *  grows beyond 2 * dx the size of the current window.           */
            w_init = tau->w_km_vals[center];

            /*  Similarly, reset the number of points in the window.          */
            nw_pts = TMPL_CAST(w_init / two_dx, size_t);
        }

        /*  Compute the Newton tranform about the current point.              */
        newton_transform(tau, x_arr, w_func, nw_pts, center);

        /*  Move to the next point. We are computing left-to-right, increment.*/
        center += 1;
    }

    /*  Free variables allocated by malloc.                                   */
    TMPL_FREE(w_func);
    TMPL_FREE(x_arr);
}
