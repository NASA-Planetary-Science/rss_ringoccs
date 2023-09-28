#include <stdlib.h>
#include <math.h>
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_special_functions_real.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      rssringoccs_Diffraction_Correction_Legendre                           *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using Legendre polynomials to           *
 *      approximate the fresnel kernel.                                       *
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
 *      1.) This routine allows for any selection of polynomial of degree     *
 *          greater than or equal to 2, though for degree 2 it is better to   *
 *          use the Fresnel option since the symmetry nearly doubles the      *
 *          speed of the computation. For anything higher than degree 8 there *
 *          is no real change in the accuracy, even for low inclination       *
 *          occultation observations.                                         *
 *      2.) Like the Fresnel approximation, the Legendre approximation has    *
 *          issues reconstructing data at low B angles. This is because the   *
 *          Legendre approximation assumes the first iteration of the Newton  *
 *          Raphson method is good enough, whereas in reality 3-4 iterations  *
 *          may be needed, like in Rev133.                                    *
 ******************************************************************************/
void rssringoccs_Diffraction_Correction_Legendre(rssringoccs_TAUObj *tau)
{
    /*  i is used for indexing, nw_pts is the number of points in the window. */
    size_t i, nw_pts, center;

    /*  Variable for the number of Legendre coefficients to be computed.      */
    unsigned int poly_order;

    /*  IsEven is a boolean for determining the parity of the polynomial.     */
    tmpl_Bool IsEven;

    /*  Various other variables needed throughout.                            */
    double w_init, dx, two_dx, cosb, sinp, cosp, Legendre_Coeff;
    double *x_arr, *w_func, *legendre_p, *alt_legendre_p, *fresnel_ker_coeffs;

    /*  Create function pointers for window function and Fresnel transform.   */
    rssringoccs_Window_Function fw = tau->window_func;

    void (*FresT)(
        rssringoccs_TAUObj *, const double *, const double *, const double *,
        size_t, size_t
    );

    /*  This should remain at false.                                          */
    tau->error_occurred = tmpl_False;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Data(tau);
    if (tau->error_occurred)
        return;

    /*  Set the IsEven boolean to the appropriate value. Since the linear and *
     *  constant term are zero, even polynomials will have an odd number of   *
     *  terms. For example, the quartic expansion has the quadratic, cubic,   *
     *  quartic terms, and hence needs three coefficients. Thus, if the       *
     *  tau->order variable is an odd number, this corresponds to an even     *
     *  polynomial and vice versa. Set IsEven accordingly.                    */
    if (tau->order & 1)
        IsEven = tmpl_True;
    else
        IsEven = tmpl_False;

    /*  Select the appropriate Fresnel transform and set poly_order.          */
    if (IsEven)
        poly_order = tau->order;
    else
        poly_order = tau->order + 1;

    if (tau->use_norm)
    {
        if (IsEven)
            FresT = rssringoccs_Fresnel_Transform_Legendre_Even_Norm;
        else
            FresT = rssringoccs_Fresnel_Transform_Legendre_Odd_Norm;
    }
    else
    {
        if (IsEven)
            FresT = rssringoccs_Fresnel_Transform_Legendre_Even;
        else
            FresT = rssringoccs_Fresnel_Transform_Legendre_Odd;
    }

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

    /*  Compute more necessary data.                                          */
    w_init = tau->w_km_vals[center];
    dx     = tau->rho_km_vals[center+1] - tau->rho_km_vals[center];
    two_dx = 2.0*dx;
    nw_pts = (size_t)(w_init / two_dx) + 1UL;

    /* Check to ensure you have enough data to the left.                      */
    rssringoccs_Tau_Check_Data_Range(tau);
    if (tau->error_occurred)
        return;

    /*  Allocate memory for the independent variable and window function.     */
    x_arr  = malloc(sizeof(*x_arr)*nw_pts);
    w_func = malloc(sizeof(*w_func)*nw_pts);

    /*  Also for the two Legendre polynomials.                                */
    legendre_p     = malloc(sizeof(*legendre_p)*(poly_order+1));
    alt_legendre_p = malloc(sizeof(*alt_legendre_p)*poly_order);

    /*  And finally for the coefficients of psi.                              */
    fresnel_ker_coeffs = malloc(sizeof(*fresnel_ker_coeffs)*poly_order);

    /*  Check that malloc was successfull then pass the x_arr array           *
     *  (ring radius) to the void function reset_window. This alters x_arr so *
     *  that it's values range from -W/2 to zero, W being the window width.   */
    if (!(x_arr)    ||    !(w_func)            ||    !(legendre_p)
                    ||    !(alt_legendre_p)    ||    !(fresnel_ker_coeffs))
    {
        /*  Malloc failed, return to calling function.                        */
        tau->error_occurred = tmpl_True;
        return;
    }
    else
        rssringoccs_Tau_Reset_Window(x_arr, w_func, dx, w_init, nw_pts, fw);

    /* Loop through each point and begin the reconstruction.                  */
    for (i = 0; i < tau->n_used; ++i)
    {
        /*  Compute some geometric information, and the scaling coefficient   *
         *  for the Legendre polynomial expansion.                            */
        cosb = tmpl_Double_Cosd(tau->B_deg_vals[center]);
        tmpl_Double_SinCosd(tau->phi_deg_vals[center], &sinp, &cosp);
        Legendre_Coeff  = cosb*sinp;
        Legendre_Coeff *= Legendre_Coeff;
        Legendre_Coeff  = 0.5*Legendre_Coeff/(1.0-Legendre_Coeff);

        /* Compute Legendre Polynomials,                                      */
        tmpl_Legendre_Polynomials(legendre_p, cosb*cosp, poly_order+1);
        tmpl_Alt_Legendre_Polynomials(alt_legendre_p, legendre_p, poly_order);

        /*  Compute the coefficients using Cauchy Products. First compute     *
         *  the bottom triangle of the square in the product.                 */
        tmpl_Fresnel_Kernel_Coefficients(fresnel_ker_coeffs, legendre_p,
                                         alt_legendre_p, Legendre_Coeff,
                                         poly_order);

        /*  If the window width changes significantly, recompute w_func.      */
        if (fabs(w_init - tau->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init = tau->w_km_vals[center];
            nw_pts = (size_t)(w_init / two_dx) + 1UL;

            /*  Reallocate x_arr and w_func since the sizes changed.          */
            x_arr  = realloc(x_arr, sizeof(double)*nw_pts);
            w_func = realloc(w_func, sizeof(double)*nw_pts);

            /*  Recompute x_arr and w_func for the new sizes.                 */
            rssringoccs_Tau_Reset_Window(x_arr, w_func, dx, w_init, nw_pts, fw);
        }

        /*  Compute the fresnel tranform about the current point.             */
        FresT(tau, x_arr, w_func, fresnel_ker_coeffs, nw_pts, center);

        /*  Increment T_in pointer using pointer arithmetic.                  */
        center += 1;
    }

    /*  Free all variables allocated by malloc.                               */
    free(x_arr);
    free(w_func);
    free(legendre_p);
    free(alt_legendre_p);
    free(fresnel_ker_coeffs);
}
