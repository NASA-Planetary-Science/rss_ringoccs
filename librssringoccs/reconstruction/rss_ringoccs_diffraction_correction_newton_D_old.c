#include <stdlib.h>
#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_string.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionNewton                                           *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the Newton-Raphson method to      *
 *      compute the stationary value of the Fresnel-Kernel.                   *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          _diffraction_correction.h. This contains all of the necessary     *
 *          data for diffraction correction, including the geometry of the    *
 *          occultation and actual power and phase data.                      *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction correction data.                  *
 *  Notes:                                                                    *
 *      1.) This method is the most accurate, but much slower than the        *
 *          Fresnel and Legendre options. It is accurate for every Rev of the *
 *          Cassini mission with the exception of the Rev133 occultation      *
 *          of which only the Ka band produces accurate results. For X and S  *
 *          bands one needs to use the Perturbed Newton method.               *
 *      2.) The polynomials from the MTR86 are available via the dlp.interp   *
 *          variable. These polynomials are slower and less accurate than the *
 *          normal Newton method since to compute the polynomials the         *
 *          Newton-Raphson method must be performed, and hence the            *
 *          polynomials increase the number of computations needed. The real  *
 *          use of them arises if one uses FFT methods. This routine does NOT *
 *          use FFTs, but rather ordinary integration.                        *
 ******************************************************************************/
void rssringoccs_Diffraction_Correction_Newton_D_Old(rssringoccs_TAUObj *tau)
{
    /*  Variables for indexing. nw_pts is the number of points in the window. */
    unsigned long i, j, nw_pts, center;

    /*  Some variables needed for reconstruction.                             */
    double w_init, dx, two_dx;
    double *x_arr, *phi_arr, *w_func;

    /*  EPS is the maximum allowed error in the Newton-Raphson scheme.        */
    double EPS = 1.E-4;

    /*  Toler is the number of iterations allowed in Newton-Raphson.          */
    unsigned long toler = 5;

    /*  Declare a function pointer for the window function.                   */
    rssringoccs_window_func fw;

    /*  If everything executes smoothly, status should remain at zero.        */
    tau->error_occurred = rssringoccs_False;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Data(tau);
    if (tau->error_occurred)
        return;

    /*  Grab the function pointer for the window function.                    */
    fw = tau->window_func;

    /* Compute first window width and window function. */
    center = tau->start;

    /*  If forward tranform is set, negate the k_vals variable. This has      *
     *  the equivalent effect of computing the forward calculation later.     */
    if (tau->use_fwd)
    {
        /*  Loop over all of k_vals and negate the value.                     */
        for (i=0; i <= tau->n_used; ++i)
            tau->k_vals[i] *= -1.0;
    }

    /*  Compute some more variables.                                          */
    w_init  = tau->w_km_vals[center];
    dx      = tau->rho_km_vals[center+1] - tau->rho_km_vals[center];
    two_dx  = 2.0*dx;
    nw_pts  = 2*((long)(w_init / two_dx))+1;

    /* Check to ensure you have enough data to the left.                      */
    rssringoccs_Tau_Check_Data_Range(tau);
    if (tau->error_occurred)
        return;

    /*  Allocate memory for these required variables.                         */
    x_arr   = malloc(sizeof(*x_arr)   * nw_pts);
    phi_arr = malloc(sizeof(*phi_arr) * nw_pts);
    w_func  = malloc(sizeof(*w_func)  * nw_pts);

    /*  Check that malloc was successfull.                                    */
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

    if (!(phi_arr))
    {
        tau->error_occurred = rssringoccs_True;
        tau->error_message = rssringoccs_strdup(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tDiffractionCorrectionFresnel\n\n"
            "\rMalloc failed and returned NULL for phi_arr. Returning.\n\n"
        );
        return;
    }

    /*  Compute the rho and phi variables, and the window function.           */
    for (j=0; j<nw_pts; ++j)
    {
        x_arr[j]   = tau->rho_km_vals[center+j-(nw_pts-1)/2];
        phi_arr[j] = tau->phi_rad_vals[center+j-(nw_pts-1)/2];
        w_func[j]  = fw(x_arr[j] - tau->rho_km_vals[center], w_init);
    }

    /*  Check tau->use_norm outside of the inner for loop to prevent          *
     *  redundantly checking an if-then statement over and over again.        */
    if (tau->use_norm)
    {
        /*  Run diffraction correction point by point.                        */
        for (i=0; i<=tau->n_used; ++i)
        {
            /*  If the window width changes significantly, recompute w_func.  */
            if (fabs(w_init - tau->w_km_vals[center]) >= two_dx)
            {
                /* Reset w_init and recompute window function.                */
                w_init  = tau->w_km_vals[center];
                nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;

                /*  Reallocate memory since the sizes have changed.           */
                w_func  = realloc(w_func,  sizeof(double) * nw_pts);
                phi_arr = realloc(phi_arr, sizeof(double) * nw_pts);
                x_arr   = realloc(x_arr,   sizeof(double) * nw_pts);

                /*  Recompute rho, phi, and the window function.              */
                for (j=0; j<nw_pts; ++j)
                {
                    x_arr[j] = tau->rho_km_vals[center+j-(nw_pts-1)/2];
                    phi_arr[j] = tau->phi_rad_vals[center+j-(nw_pts-1)/2];
                    w_func[j] = fw(x_arr[j] - tau->rho_km_vals[center], w_init);
                }
            }
            else
            {
                /*  Adjust rho and phi to the new range.                      */
                for (j=0; j<nw_pts; ++j)
                {
                    x_arr[j]   = tau->rho_km_vals[center+j-(nw_pts-1)/2];
                    phi_arr[j] = tau->phi_rad_vals[center+j-(nw_pts-1)/2];
                }
            }

            /*  Compute the fresnel tranform about the current point.         */
            tau->T_out[i] = Fresnel_Transform_Newton_D_Old_Norm_Double(
                x_arr, phi_arr, tau->T_in, w_func, tau->k_vals[center],
                tau->rho_km_vals[center], tau->B_rad_vals[center],EPS, toler,
                nw_pts, center, tau->rx_km_vals[center],
                tau->ry_km_vals[center], tau->rz_km_vals[center]
            );

            /*  Increment pointers using pointer arithmetic.                  */
            center += 1;
        }
    }
    else
    {
        /*  Run diffraction correction point by point.                        */
        for (i=0; i<=tau->n_used; ++i)
        {
            /*  If the window width changes significantly, recompute w_func.  */
            if (fabs(w_init - tau->w_km_vals[center]) >= two_dx)
            {
                /* Reset w_init and recompute window function.                */
                w_init  = tau->w_km_vals[center];
                nw_pts  = 2*((int)(w_init / (2.0 * dx)))+1;

                /*  Reallocate memory since the sizes have changed.           */
                w_func  = realloc(w_func,  sizeof(double) * nw_pts);
                phi_arr = realloc(phi_arr, sizeof(double) * nw_pts);
                x_arr   = realloc(x_arr,   sizeof(double) * nw_pts);

                /*  Recompute rho, phi, and the window function.              */
                for (j=0; j<nw_pts; ++j)
                {
                    x_arr[j] = tau->rho_km_vals[center+j-(nw_pts-1)/2];
                    phi_arr[j] = tau->phi_rad_vals[center+j-(nw_pts-1)/2];
                    w_func[j] = fw(x_arr[j] - tau->rho_km_vals[center], w_init);
                }
            }
            else
            {
                /*  Adjust rho and phi to the new range.                      */
                for (j=0; j<nw_pts; ++j)
                {
                    x_arr[j]   = tau->rho_km_vals[center+j-(nw_pts-1)/2];
                    phi_arr[j] = tau->phi_rad_vals[center+j-(nw_pts-1)/2];
                }
            }

            /*  Compute the fresnel tranform about the current point.         */
            tau->T_out[i] = Fresnel_Transform_Newton_D_Old_Double(
                x_arr, phi_arr, tau->T_in, w_func, tau->k_vals[center],
                tau->rho_km_vals[center], tau->B_rad_vals[center],
                EPS, toler, dx, tau->F_km_vals[center], nw_pts, center,
                tau->rx_km_vals[center], tau->ry_km_vals[center],
                tau->rz_km_vals[center]
            );

            /*  Increment pointers using pointer arithmetic.                  */
            center += 1;
        }
    }

    /*  Free variables allocated by malloc.                                   */
    free(x_arr);
    free(phi_arr);
    free(w_func);
}
