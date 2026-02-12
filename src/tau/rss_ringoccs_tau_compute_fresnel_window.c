#include <libtmpl/include/tmpl_config.h>
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/constants/tmpl_math_constants.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

void
rssringoccs_Tau_Reset_Fresnel_Window(
    const rssringoccs_TAUObj * TMPL_RESTRICT const tau,
    double * TMPL_RESTRICT const x_arr,
    double * TMPL_RESTRICT const w_func,
    const size_t nw_pts,
    const size_t center
)
{
    /*  Create a variable for indexing.                                       */
    size_t n;

    /*  The nth element in the x array is (n - nw) * dx, cast to double. We   *
     *  perform the cast first since size_t is unsigned, hence if n - nw is   *
     *  negative, the computation will wrap around to a very large number and *
     *  give us a gibberish result.                                           */
    const double offset = TMPL_CAST(nw_pts, double);

    /* Loop over n, computing the window function and the x_arr variable.     */
    for (n = 0; n < nw_pts; ++n)
    {
        const double index = TMPL_CAST(n, double);
        x_arr[n] = (index - offset) * tau->dx_km;
        w_func[n] = tau->window_func(x_arr[n], tau->w_km_vals[center]);

        /*  We have computed the window function and the independent variable *
         *  x, which is (r - r0). We need (pi/2) (r - r0)^2. The 1 / F^2      *
         *  factor is introduced later inside the Fresnel transform function. */
        x_arr[n] *= tmpl_double_pi_by_two * x_arr[n];
    }
}
