#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

int
rssringoccs_Tau_Resize_Half_Window(
    rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    double ** TMPL_RESTRICT const x_ptr,
    double ** TMPL_RESTRICT const w_ptr,
    const double width,
    const double two_dx,
    const size_t center
)
{
    /*  Do not attempt to access a NULL pointer. Check for this.              */
    if (!tau)
        return 1;

    /*  If an error occurred before calling this function, abort.             */
    if (tau->error_occurred)
        return 1;

    /*  More checks for NULL pointers.                                        */
    if (!x_ptr || !w_ptr)
    {
        /*  Exit with an error message.                                       */
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Tau_Check_And_Resize_Window\n\n"
            "\rInput pointers are NULL. Aborting.\n\n";

        return 1;
    }

    /*  If the window width changes significantly, recompute window array.    */
    if (tmpl_Double_Abs(width - tau->w_km_vals[center]) >= two_dx)
    {
        /*  Temporary pointers for realloc. This will help avoid leaks.       */
        void *x_tmp, *w_tmp;

        /*  Pointers to the arrays for x and the window function.             */
        double *x_arr = *x_ptr;
        double *w_func = *w_ptr;

        /*  Variables for the window size (in km) and half the number of      *
         *  points that fit inside the window.                                */
        double w_init;
        size_t nw_pts;

        /*  Check that the arrays that x_ptr and w_ptr are pointing to are    *
         *  not NULL. Most routines that call this function allocate memory   *
         *  for these arrays via malloc. It is possible malloc failed and     *
         *  this function recieved pointers to NULL pointers. Check for this. */
        if (!x_arr || !w_func)
        {
            /*  Exit with an error message.                                   */
            tau->error_occurred = tmpl_True;
            tau->error_message = tmpl_String_Duplicate(
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Check_And_Resize_Window\n\n"
                "\rInput pointers are pointing to NULL pointers. Aborting.\n\n"
            );

            return 1;
        }

        /*  Reset the window size and recompute the number of points needed.  */
        w_init = tau->w_km_vals[center];
        nw_pts = TMPL_CAST(w_init / two_dx, size_t);

        /*  Reallocate memory for the x array since the sizes changed.        */
        x_tmp = realloc(x_arr, sizeof(*x_arr) * nw_pts);

        /*  Check to make sure realloc didn't fail.                           */
        if (!x_tmp)
        {
            /*  Abort the computation with an error message.                  */
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Check_And_Resize_Window\n\n"
                "\rrealloc failed and returned NULL for x_arr.\n\n";

            /*  If we get here, realloc failed for x_arr. This means we can   *
             *  safely free x_arr. x_tmp is just a NULL pointer, do not       *
             *  attempt to free it.                                           */
            TMPL_FREE(x_arr);

            /*  We have not attempted to realloc w_func yet, meaning this     *
             *  pointer still points to the initial data. Free this.          */
            TMPL_FREE(w_func);

            /*  The memory for x_arr and w_func has been freed, set the       *
             *  input pointers to NULL to avoid double free's later.          */
            *x_ptr = NULL;
            *w_ptr = NULL;

            /*  realloc failed, which is treated as an error. Return 1.       */
            return 1;
        }

        /*  Similarly, reallocate memory for the window function.             */
        w_tmp = realloc(w_func, sizeof(*w_func) * nw_pts);

        /*  Check to make sure realloc didn't fail.                           */
        if (!w_tmp)
        {
            /*  Abort the computation with an error message.                  */
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Check_And_Resize_Window\n\n"
                "\rrealloc failed and returned NULL for w_func.\n\n";

            /*  realloc succeeded with x_arr, but failed with w_func. The     *
             *  temporary pointer hence points to the valid data, free this.  */
            TMPL_FREE(x_tmp);

            /*  realloc failed for w_tmp, meaning w_tmp is NULL and w_func    *
             *  points to valid data. Free this.                              */
            TMPL_FREE(w_func);

            /*  The memory for x_arr and w_func has been freed, set the       *
             *  input pointers to NULL to avoid double free's later.          */
            *x_ptr = NULL;
            *w_ptr = NULL;

            /*  realloc failed, which is treated as an error. Return 1.       */
            return 1;
        }

        /*  If we get here, realloc succeeded. Swap the variables.            */
        *x_ptr = x_tmp;
        *w_ptr = w_tmp;

        /*  Recompute x_ptr and w_ptr for the new sizes.                      */
        rssringoccs_Tau_Reset_Window(
            tau, *x_ptr, *w_ptr, nw_pts, center
        );

        /*  Window resize was needed, return 1.                               */
        return 1;
    }

    /*  No error occurred and there is no need to resize the window. Return 0.*/
    return 0;
}
