#include <libtmpl/include/tmpl_config.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

void
rssringoccs_Tau_Compute_Window(rssringoccs_TAUObj * TMPL_RESTRICT const tau,
                               double * TMPL_RESTRICT const w_func,
                               const size_t nw_pts,
                               const size_t center)
{
    /*  The index for the offset (the left-most edge of the window) is just   *
     *  half the number of points in the window away to the left. Compute.    */
    const size_t offset = center - (nw_pts >> 1);

    /*  Create a variable for indexing.                                       */
    size_t n;

    /*  Compute the values in the window function.                            */
    for (n = 0; n < nw_pts; ++n)
    {
        /*  We compute left-to-right. The offset has the index of the left    *
         *  most point, we can use this to obtain the index for the current   *
         *  point for the window function.                                    */
        const size_t ind = offset + n;

        /*  The tapering function is a function of (rho - rho0), where rho is *
         *  the radius of the current point, and rho0 is the radius of the    *
         *  central point. Compute this.                                      */
        const double x = tau->rho_km_vals[ind] - tau->rho_km_vals[center];

        /*  The Tau object contains a function pointer to the selected window *
         *  function. Compute using this.                                     */
        w_func[n] = tau->window_func(x, tau->w_km_vals[center]);
    }
}
