#include <libtmpl/include/compat/tmpl_cast.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      reset_window                                                          *
 *  Purpose:                                                                  *
 *      Compute an array of points from -width to zero, equally spaced by dx. *
 *      This acts as the independent variable for later use.                  *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed, and rho is the dummy variable of integration *
 *          which varies from rho0-W/2 to rho0, W being the window width.     *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr.             *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      width (double):                                                       *
 *          The width of the window function.                                 *
 *      nw_pts (long):                                                        *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *          with only half of the window. This saves a lot of computation.    *
 *      fw (rss_ringoccs_window_func):                                        *
 *          Function pointer to the window function.                          *
 *  Notes:                                                                    *
 *      1.) This is a void function that takes in pointers as arguments. The  *
 *          values of the pointers are changed within this function and there *
 *          is no need to return anything. Hence, no return statement.        *
 *      2.) x_arr, width, and dx should all be in the same units. It is       *
 *          recommended one uses kilometers. Constants like the speed of      *
 *          are defined using km/s, and the CSV files which contain the       *
 *          occultation data store values in kilometers.                      *
 ******************************************************************************/
void
rssringoccs_Tau_Reset_Window(double * TMPL_RESTRICT const x_arr,
                             double * TMPL_RESTRICT const w_func,
                             double dx,
                             double width,
                             size_t nw_pts,
                             rssringoccs_WindowFunction window)
{
    /*  Create a variable for indexing.                                       */
    size_t n;

    /*  The nth element in the x array is (n - nw) * dx, cast to double. We   *
     *  perform the cast first since size_t is unsigned, hence if n - nw is   *
     *  negative, the computation will wrap around to a very large number and *
     *  give us a gibberish result.                                           */
    const double offset = TMPL_CAST(nw_pts, double);

    /* Loop over n, computing the window function and the x_arr variable.     */
    for(n = 0U; n < nw_pts; ++n)
    {
        const double index = TMPL_CAST(n, double);
        x_arr[n] = (index - offset)*dx;
        w_func[n] = window(x_arr[n], width);
    }
}
