#include <stdlib.h>
#include <math.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_string.h>
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
void rssringoccs_Diffraction_Correction_Newton(rssringoccs_TAUObj *tau)
{
    /*  Variables for indexing. nw_pts is the number of points in the window. */
    size_t i, j, offset, nw_pts, center;

    /*  Some variables needed for reconstruction.                             */
    double w_init, dx, two_dx;
    double *w_func;

    /*  Declare a function pointer for the transform function.                */
    rssringoccs_FresT FresT;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Data(tau);

    if (tau->error_occurred)
        return;

    /* Compute first window width and window function. */
    center = tau->start;

    /*  Compute some more variables.                                          */
    w_init = tau->w_km_vals[center];
    dx     = tau->rho_km_vals[center+1] - tau->rho_km_vals[center];
    two_dx = 2.0*dx;
    nw_pts = 2UL*((size_t)(w_init / two_dx)) + 1UL;
    offset = center - (nw_pts - 1UL) / 2UL;

    /* Check to ensure you have enough data to the left.                      */
    rssringoccs_Tau_Check_Data_Range(tau);

    if (tau->error_occurred)
        return;

    /*  Set the correct function pointer.                                     */
    if (tau->use_norm)
    {
        if (tau->psinum == rssringoccs_DR_Newton)
            FresT = rssringoccs_Fresnel_Transform_Newton_Norm;
        else if (tau->psinum == rssringoccs_DR_NewtonD)
            FresT = rssringoccs_Fresnel_Transform_Newton_D_Norm;
        else if (tau->psinum == rssringoccs_DR_NewtonDPhi)
            FresT = rssringoccs_Fresnel_Transform_Newton_dD_dphi_Norm;
        else if (tau->psinum == rssringoccs_DR_NewtonPerturb)
            FresT = rssringoccs_Fresnel_Transform_Perturbed_Newton_Norm;
        else if (tau->psinum == rssringoccs_DR_NewtonElliptical)
            FresT = rssringoccs_Fresnel_Transform_Newton_Elliptical_Norm;
        else if (tau->psinum == rssringoccs_DR_NewtonQuadratic)
            FresT = rssringoccs_Fresnel_Transform_Newton_Quadratic_Norm;
        else if (tau->psinum == rssringoccs_DR_NewtonQuartic)
            FresT = rssringoccs_Fresnel_Transform_Newton_Quartic_Norm;
        else if (tau->psinum == rssringoccs_DR_NewtonD)
            FresT = rssringoccs_Fresnel_Transform_Newton_D_Quartic_Norm;
        else
            FresT = rssringoccs_Fresnel_Transform_Newton_D_Old_Norm;
    }
    else
    {
        if (tau->psinum == rssringoccs_DR_Newton)
            FresT = rssringoccs_Fresnel_Transform_Newton;
        else if (tau->psinum == rssringoccs_DR_NewtonD)
            FresT = rssringoccs_Fresnel_Transform_Newton_D;
        else if (tau->psinum == rssringoccs_DR_NewtonDPhi)
            FresT = rssringoccs_Fresnel_Transform_Newton_dD_dphi;
        else if (tau->psinum == rssringoccs_DR_NewtonPerturb)
            FresT = rssringoccs_Fresnel_Transform_Perturbed_Newton;
        else if (tau->psinum == rssringoccs_DR_NewtonElliptical)
            FresT = rssringoccs_Fresnel_Transform_Newton_Elliptical;
        else if (tau->psinum == rssringoccs_DR_NewtonQuadratic)
            FresT = rssringoccs_Fresnel_Transform_Newton_Quadratic;
        else if (tau->psinum == rssringoccs_DR_NewtonQuartic)
            FresT = rssringoccs_Fresnel_Transform_Newton_Quartic;
        else if (tau->psinum == rssringoccs_DR_NewtonD)
            FresT = rssringoccs_Fresnel_Transform_Newton_D_Quartic;
        else
            FresT = rssringoccs_Fresnel_Transform_Newton_D_Old;
    }

    /*  Allocate memory for these required variables.                         */
    w_func = malloc(sizeof(*w_func) * nw_pts);

    if (!(w_func))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tDiffractionCorrectionFresnel\n\n"
            "\rMalloc failed and returned NULL for w_func. Returning.\n\n"
        );
        return;
    }

    /*  Compute the rho and phi variables, and the window function.           */
    for (j=0; j<nw_pts; ++j)
        w_func[j] = tau->window_func(
            tau->rho_km_vals[offset+j] - tau->rho_km_vals[center], w_init
        );

    /*  Check tau->use_norm outside of the inner for loop to prevent          *
     *  redundantly checking an if-then statement over and over again.        */

    /*  Run diffraction correction point by point.                        */
    for (i=0; i < tau->n_used; ++i)
    {
        /*  If the window width changes significantly, recompute w_func.  */
        if (fabs(w_init - tau->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                */
            w_init = tau->w_km_vals[center];
            nw_pts = 2*((size_t)(w_init / (2.0 * dx))) + 1UL;
            offset = center - (nw_pts - 1UL) / 2UL;

            /*  Reallocate memory since the sizes have changed.           */
            w_func = realloc(w_func, sizeof(*w_func) * nw_pts);

            /*  Recompute rho, phi, and the window function.              */
            for (j=0; j<nw_pts; ++j)
                w_func[j] = tau->window_func(
                    tau->rho_km_vals[offset+j] - tau->rho_km_vals[center],
                    w_init
                );
        }

        /*  Compute the fresnel tranform about the current point.         */
        FresT(tau, w_func, nw_pts, center);

        /*  Increment pointers using pointer arithmetic.                  */
        center += 1;
    }

    /*  Free variables allocated by malloc.                                   */
    free(w_func);
}
