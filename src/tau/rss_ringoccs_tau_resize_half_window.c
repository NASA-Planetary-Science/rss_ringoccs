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
        double * x_arr = *x_ptr;
        double * w_func = *w_ptr;

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

        /*  Reallocate memory for x and w since the sizes changed.            */
        x_tmp = realloc(x_arr, sizeof(*x_arr) * nw_pts);
        w_tmp = realloc(w_func, sizeof(*w_func) * nw_pts);

        /*  Check to make sure realloc didn't fail.                           */
        if (!x_tmp || !w_tmp)
        {
            /*  Abort the computation with an error message.                  */
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Tau_Check_And_Resize_Window\n\n"
                "\rrealloc failed and returned NULL. Aborting.\n\n";

            /*  Free all variables allocated by malloc. Since we only know    *
             *  realloc failed for either x_arr or w_func, but not which one, *
             *  we need to check before freeing. Freeing a pointer after a    *
             *  successful call to realloc is undefined behavior. First check *
             *  if x_tmp passed. If it did, free x_tmp instead of x_ptr.      */
            TMPL_FREE(x_tmp);

            /*  If we get here, realloc failed for x_arr. This means we can   *
             *  safely free x_arr. x_tmp is just a NULL pointer, do not       *
             *  attempt to free it.                                           */
            TMPL_FREE(x_arr);

            /*  Same check for w_func. It is possible that realloc failed for *
             *  x_ptr, but worked for w_func. Checking w_tmp tells us the     *
             *  answer. If w_tmp is NULL, realloc failed for w_func meaning   *
             *  we should free w_func. If w_tmp is not NULL, then realloc was *
             *  successful and we must free w_tmp. Check this.                */
            TMPL_FREE(w_tmp);

            /*  If we get here, realloc failed for w_func. We may safely free *
             *  w_func. w_tmp is just a NULL pointer, do not free it.         */
            TMPL_FREE(w_func);

            /*  At this point, the memory allocated to w_ptr and x_ptr have   *
             *  been freed. Set these to NULL to prevent the possibility of   *
             *  freeing these pointers twice.                                 */
            *x_ptr = NULL;
            *w_ptr = NULL;

            return 1;
        }

        /*  If we get here, realloc succeeded. Swap the variables.            */
        *x_ptr = x_tmp;
        *w_ptr = w_tmp;

        /*  Recompute x_ptr and w_ptr for the new sizes.                      */
        rssringoccs_Tau_Reset_Window(
            tau, *x_ptr, *w_ptr, nw_pts, center
        );

        return 1;
    }

    /*  No error occurred and there is no need to resize the window. Return 0.*/
    return 0;
}
