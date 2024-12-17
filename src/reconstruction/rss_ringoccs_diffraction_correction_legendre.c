/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  rss_ringoccs is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************
 *              rss_ringoccs_diffraction_correction_legendre_fast             *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Performs diffraction correction using Legendre polynomials, and by    *
 *      assuming rho / rho_0 ~= 1 across the window of integration.           *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Diffraction_Correction_Legendre_Fast                      *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using Legendre polynomials, and using *
 *      the approximation rho / rho_0 ~= 1.                                   *
 *  Arguments:                                                                *
 *      tau (double * const):                                                 *
 *          The geometry and diffraction data for the reconstruction. The     *
 *          output reconstruction is stored in the T_out member of tau.       *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      None.                                                                 *
 *  Method:                                                                   *
 *      The Fresnel transform is given by:                                    *
 *  Notes:                                                                    *
 *      1.) This function assumes tau->order > 1. Degree zero and degree one  *
 *          reconstructions are not allowed since these polynomials are zero. *
 *          Degree 2 is valid, but this is the classic Fresnel approximation, *
 *          and one should use rssringoccs_Diffraction_Correction_Fresnel.    *
 *          For degree greater than 2, we use the expansion found in [1].     *
 *      2.) This function is not accurate for very large windows since the    *
 *          assumption rho / rho_0 ~= 1 is no longer valid. Be wary of this.  *
 *      3.) For low opening angles this function yields a poor reconstruction.*
 *          You are better off using the Newton-Raphson based routines.       *
 *  References:                                                               *
 *      1.) Maguire, R., French, R. (2024)                                    *
 *          "Applications of Legendre Polynomials for Fresnel Inversion       *
 *              and Occultation Observations"                                 *
 *      2.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          "Profiling Saturn's Rings by Radio Occultation"                   *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *      3.) Goodman, J. (2005)                                                *
 *          "Introduction to Fourier Optics"                                  *
 *          McGraw-Hill Series in Electrical and Computer Engineering.        *
 *      4.) McQuarrie, Donald (2003),                                         *
 *          "Mathematical Methods for Scientists and Engineers",              *
 *          University Science Books, ISBN 1-891389-29-7,                     *
 *          Chapter 14 "Orthogonal Polynomials and Sturm-Liouville Problems"  *
 *      5.) Arfken, G., Weber, H., Harris, F. (2013)                          *
 *          "Mathematical Methods for Physicists, Seventh Edition"            *
 *          Academic Press, Elsevier                                          *
 *          Chapter 15 "Legendre Functions"                                   *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *                                                                            *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       December 11, 2024                                             *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2024/12/14: Ryan Maguire                                                  *
 *      Added license, fixed comments.                                        *
 ******************************************************************************/
#include <libtmpl/include/tmpl_math.h>
#include <libtmpl/include/tmpl_orthogonal_polynomial_real.h>
#include <libtmpl/include/tmpl_compat_cast.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stdlib.h>

static const rssringoccs_FresnelLegendreTransform
legendre_transform_list[4] = {
    rssringoccs_Fresnel_Transform_Legendre_Odd,
    rssringoccs_Fresnel_Transform_Legendre_Even,
    rssringoccs_Fresnel_Transform_Legendre_Odd_Norm,
    rssringoccs_Fresnel_Transform_Legendre_Even_Norm
};

#define RSSRINGOCCS_DESTROY_VARIABLE(var) if (var) free(var)

void rssringoccs_Diffraction_Correction_Legendre(rssringoccs_TAUObj *tau)
{
    /*  n is used for indexing, nw_pts is the number of points in the window. */
    size_t n, nw_pts, center, index;

    /*  is_even is a boolean for determining the parity of the polynomial.    */
    tmpl_Bool is_even;

    /*  Various other variables needed throughout.                            */
    double w_init, dx, two_dx, cosb, sinp, cosp, alpha, beta;

    /*  Pointers for arrays. Initialize them to NULL so that we may easily    *
     *  check if malloc fails and avoid freeing non-malloced pointers.        */
    double *x_arr = NULL;
    double *w_func = NULL;
    double *coeffs = NULL;

    /*  Create function pointers for window function and Legendre transforms. */
    rssringoccs_WindowFunction fw;
    rssringoccs_FresnelLegendreTransform legendre_transform;

    /*  Make sure the input is not NULL before checking its data.             */
    if (!tau)
        return;

    /*  If an error occurred before we got to this function, abort.           */
    if (tau->error_occurred)
        return;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Data(tau);

    /* Check to ensure you have enough data to process.                       */
    rssringoccs_Tau_Check_Data_Range(tau);

    /*  The above routines set the error_occurred flag on error. Check.       */
    if (tau->error_occurred)
        return;

    /*  The "is_even" Boolean tells us whether the degree of the polynomial   *
     *  is even, not the number of coefficients. Thus, if a quartic           *
     *  Legendre polynomial is desired, is_even will be true. The number of   *
     *  coefficients needed has the opposite parity. For example, a quartic   *
     *  polynomial needs the quadratic, cubic, and quartic terms since the    *
     *  linear and constant terms are zero. Hence there are three             *
     *  coefficients, which is odd. tau->order has the number of coefficients *
     *  needed, and not the degree of the polynomial. The degree is order + 1.*
     *  We can determine if the degree is even by checking if order is odd.   *
     *  We check if order is odd by examining the last bit in tau->order. If  *
     *  it is 1, tau->order is odd (so degree is even), and if it is zero we  *
     *  have that tau->order is even (so degree is odd).                      */
    is_even = (tau->order & 0x01U);

    /*  There are four types of Legendre functions:                           *
     *      0.) Odd Degree without Window Normalization.                      *
     *      1.) Even Degree without Window Normalization.                     *
     *      2.) Odd Degree with Window Normalization.                         *
     *      3.) Even Degree with Window Normalization.                        *
     *  The use_norm member is a Boolean, zero or one, as is the is_even      *
     *  variable. The index above can be compute by is_even + 2*use_norm.     */
    index = TMPL_CAST(is_even, size_t) + 2 * TMPL_CAST(tau->use_norm, size_t);
    legendre_transform = legendre_transform_list[index];

    /*  If forward tranform is set, negate the k_vals variable. This has      *
     *  the equivalent effect of computing the forward calculation later.     */
    if (tau->use_fwd)
    {
        /*  Loop over all of k_vals and negate the value.                     */
        for (n = 0; n <= tau->n_used; ++n)
            tau->k_vals[n] *= -1.0;
    }

    /*  Compute necessary data for the start of the inversion.                */
    center = tau->start;
    w_init = tau->w_km_vals[center];
    dx = tau->rho_km_vals[center + 1] - tau->rho_km_vals[center];
    two_dx = 2.0 * dx;
    nw_pts = TMPL_CAST(w_init / two_dx, size_t) + 1;
    fw = tau->window_func;

    /*  Allocate memory for the independent variable and window function.     */
    x_arr = malloc(sizeof(*x_arr) * nw_pts);
    w_func = malloc(sizeof(*w_func) * nw_pts);

    /*  And finally for the coefficients of the Fresnel-Legendre polynomial.  */
    coeffs = malloc(sizeof(*coeffs) * tau->order);

    /*  Check if malloc failed. It returns NULL on failure.                   */
    if (!x_arr || !w_func || !coeffs)
    {
        /*  malloc failed. Return with error, and free any successfully       *
         *  allocated data.                                                   */
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_String_Duplicate(
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Diffraction_Correction_Fast_Legendre\n\n"
            "\rmalloc failed and returned NULL. Aborting.\n\n"
        );

        /*  It is possible that malloc succeeded for some variables, and not  *
         *  for others. Since we initialized these pointers to NULL, if a     *
         *  pointer is not NULL, then malloc was successful in this case and  *
         *  we need to free it. The RSSRINGOCCS_DESTROY_VARIABLE macro only   *
         *  calls free if the pointer is not NULL. Use this on each variable. */
        RSSRINGOCCS_DESTROY_VARIABLE(x_arr);
        RSSRINGOCCS_DESTROY_VARIABLE(w_func);
        RSSRINGOCCS_DESTROY_VARIABLE(coeffs);

        return;
    }

    /*  Initialize the window function and the indpendent variable "x".       */
    rssringoccs_Tau_Reset_Window(x_arr, w_func, dx, w_init, nw_pts, fw);

    /* Loop through each point and begin the reconstruction.                  */
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
            {
                /*  The x_arr and w_func data are handled by the previous     *
                 *  function call on error. We still need to free the data    *
                 *  allocated for the Fresnel coefficients.                   */
                free(coeffs);
                return;
            }

            /*  Reset the threshold value for the window width to the current *
             *  window size. We will update again once the required window    *
             *  grows beyond 2 * dx the size of the current window.           */
            w_init = tau->w_km_vals[center];

            /*  Similarly, reset the number of points in the window.          */
            nw_pts = TMPL_CAST(w_init / two_dx, size_t) + 1;
        }

        /*  Compute some geometric information, and the scaling coefficient   *
         *  for the Legendre polynomial expansion.                            */
        cosb = tmpl_Double_Cosd(tau->B_deg_vals[center]);
        tmpl_Double_SinCosd(tau->phi_deg_vals[center], &sinp, &cosp);
        alpha = cosb * cosp;
        beta = cosb * sinp;
        beta *= beta;
        beta *= 0.5 / (1.0 - beta);

        /*  Initialize T_out to zero so we can loop over later.               */
        tau->T_out[center] = tmpl_CDouble_Zero;

        /*  The Fresnel-Legendre coefficients are computed using an upwards   *
         *  recursion in terms of the Legendre polynomials, and the Chebyshev *
         *  polynomials of the second kind. We have:                          *
         *                                                                    *
         *                 P (a) - a P   (a)      -                     -     *
         *                  n         n+1        |                       |    *
         *      L (a, b) = ----------------- - b |  U   (a) - 2 P   (a)  |    *
         *       n               n + 2           |   n+2         n+2     |    *
         *                                        -                     -     *
         *                                                                    *
         *  Where Pn is the nth Legendre polynomial, Un is the nth Chebyshev  *
         *  polynomial of the second kind, and we have used alpha = a and     *
         *  beta = b in order for the equation to fit onto the screen. This   *
         *  is computed by libtmpl in upwards iterative fashion.              */
        tmpl_Double_Fresnel_Legendre_L(coeffs, alpha, beta, tau->order);

        /*  Compute the Legendre tranform about the current point.            */
        legendre_transform(tau, x_arr, w_func, coeffs, nw_pts, center);

        /*  Move on to the next point in the data.                            */
        center += 1;
    }

    /*  Free all variables allocated by malloc.                               */
    free(x_arr);
    free(w_func);
    free(coeffs);
}

/*  Undefine everything in case someone wants to #include this file.          */
#undef RSSRINGOCCS_DESTROY_VARIABLE
