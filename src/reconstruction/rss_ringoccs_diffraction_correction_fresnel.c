/******************************************************************************
 *                                 LICENSE                                    *
 ******************************************************************************
 *  This file is part of rss_ringoccs.                                        *
 *                                                                            *
 *  rss_ringoccs is free software: you can redistribute it and/or modify it   *
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
 *                          Diffraction Correction                            *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      This file contains routines to compute the Fresnel transform for      *
 *      diffraction reconstruction or forward modeling.                       *
 *          Fresnel Quadratic Approximations:                                 *
 *              Classic quadratic approximation used in Fourier optics.       *
 *          Legendre Expansion:                                               *
 *             Approximate the Fresnel kernel using Legendre polynomials.     *
 *          Newton-Raphson Method:                                            *
 *              Performs the Fresnel Transform by computing the stationary    *
 *              value of the Fresnel Kernel using the Newton-Raphson method.  *
 *          Simple FFT Method:                                                *
 *              Uses an FFT to compute the Fresnel transform. This assumes    *
 *              the geometry about the center of the data set is an accurate  *
 *              representative of the entire occultation. Hence the method is *
 *              accurate around the center and inaccurate towards the edges.  *
 ******************************************************************************
 *  The Inverse Fresnel Transform:                                            *
 *                                                                            *
 *                           infinity                                         *
 *                              -                                             *
 *                    1 + i    | |                                            *
 *         T(rho)  =  ------   |   T_hat(r_0)w(r-r_0)exp(-i psi(r,r_0)) dr_0  *
 *                      2F   | |                                              *
 *                            -                                               *
 *                          -infinity                                         *
 *                                                                            *
 *      Where T_hat is the diffracted data, w is the window function, r is    *
 *      the ring intercept point, and r_0 is a dummy variable of integration. *
 *      psi is the Fresnel Kernel, and exp is the exponential function. The   *
 *      scale factor F is the Fresnel scale which is dependent on the         *
 *      geometry of the occultation.                                          *
 *                                                                            *
 *      For near-constant geometry, or regions where the geometry does not    *
 *      depend heavily on phi, psi is approximately a function of the form    *
 *      psi(rho,rho_0) = f(rho-rho_0) and so this integral is a convolution   *
 *      of w * exp(-i psi) and T_hat. Using the convolution theorem, we get:  *
 *                                                                            *
 *      F(T)   =    F(T_hat) * F(w * exp(-i psi))                             *
 *                                                                            *
 *      T      =    F^-1(F(T_hat) * F(w * exp(-i psi)))                       *
 *                                                                            *
 *      Where F(f) is the Fourier transform of f, and F^-1 is the inverse     *
 *      Fourier transform.                                                    *
 ******************************************************************************
 *  The Normalization Scheme:                                                 *
 *      As the resolution get's too high, say 10 km or larger, the window     *
 *      width quickly shrinks to zero and the integral will be approximately  *
 *      zero. To account for this, the option to normalize the integral by    *
 *      the window width is offered. The normalization is defined as follows: *
 *                                                                            *
 *                    |     _ +infinity           |                           *
 *                    |    | |                    |                           *
 *                    |    |    exp(-i psi(x)) dx |                           *
 *                    |  | |                      |                           *
 *                    |   -  -infinity            |                           *
 *          Norm =  __________________________________                        *
 *                  |    -  +W/2                    |                         *
 *                  |   | |                         |                         *
 *                  |   |    w(x) exp(-i psi(x)) dx |                         *
 *                  | | |                           |                         *
 *                  |  -   -W/2                     |                         *
 *                                                                            *
 *      This has the effect of making the free-space regions, which are       *
 *      regions that are not affected by diffraction, evaluate to one,        *
 *      regardless of what (positive) resolution was chosen.                  *
 ******************************************************************************
 *                              DEFINED FUNCTIONS                             *
 ******************************************************************************
 *  reset_window:                                                             *
 *      Void function that takes in a pointer to a double array and creates   *
 *      an array of values for the ring radius within half a window width of  *
 *      the ring intercept point. Do to symmetry, only the values to the left *
 *      of the ring intercept point are computed. This also computes the      *
 *      window function as a function of this array.                          *
 ******************************************************************************
 *  Fresnel_Transform_Double / Fresnel_Transform_Norm_Double                  *
 *      Computes the Fresnel transform using a quadratic approximation to the *
 *      Fresnel kernel, with or without normalization.                        *
 ******************************************************************************
 *  Fresnel_Legendre_Double / Fresnel_Legendre_Norm_Double                    *
 *      Computes the Fresnel Transform using Legendre Polynomials to          *
 *      approximate the Fresnel kernel to various powers (1<n<256).           *
 ******************************************************************************
 *  Fresnel_Transform_Newton_Double / Fresnel_Transform_Newton_Norm_Double    *
 *      Computes the Fresnel inverse transform using Newton-Raphson to        *
 *      compute the stationary value of the Fresnel kernel.                   *
 ******************************************************************************
 *                               DEPENDENCIES                                 *
 ******************************************************************************
 *  1.) complex.h:                                                            *
 *      A standard library header file where complex types are defined and    *
 *      various functions for manipulating complex values. This is available  *
 *      in the C99 standard and higher (C11, C18).                            *
 *  2.) string.h:                                                             *
 *      A standard library header file used for dealing with strings. This is *
 *      primarily used for checking the psitype and wtype inputs.             *
 *  3.) stdbool.h:                                                            *
 *      A standard library header file used for dealing with Booleans. It     *
 *      provides the alias bool for the built-in _Bool type. It also provides *
 *      true and false macros. Requires C99 or higher.                        *
 ******************************************************************************
 *                            A NOTE ON COMMENTS                              *
 ******************************************************************************
 *  It is anticipated that many users of this code will have experience in    *
 *  either Python or IDL, but not C. Many comments are left to explain as     *
 *  much as possible. Vagueness or unclear code should be reported to:        *
 *  https://github.com/NASA-Planetary-Science/rss_ringoccs/issues             *
 ******************************************************************************
 *                          A NOTE ON CONVENTIONS                             *
 ******************************************************************************
 *  1.) i is a complex number, we'll use n for indexing. Max index values for *
 *      a for loop or while loop should be either capital N, or preferrably   *
 *      a variable name describing what this limit is. i.e. if it's the size  *
 *      of an array, call it arr_size.                                        *
 *  2.) Do not cast malloc. While this is required in C++, it is not in C and *
 *      the official documentation even frowns upon it. That is, do this:     *
 *          double *ptr;                                                      *
 *          long N;                                                           *
 *          ...                                                               *
 *          ptr = malloc(sizeof(*ptr) * N);                                   *
 *      and not:                                                              *
 *          ptr = (double *)malloc(sizeof(*ptr) * N);                         *
 *      malloc returns a void pointer and this is automatically and safely    *
 *      promoted to whatever the type of ptr is.                              *
 *  3.) While not frowned upon, the following makes maintenance easier.       *
 *      Instead of using malloc like (for a pointer to a double):             *
 *          ptr = malloc(sizeof(double) * N);                                 *
 *      Use:                                                                  *
 *          ptr = malloc(sizeof(*ptr) * N);                                   *
 *  4.) Declare variables towards the top of the code and not inside a for    *
 *      loop. It's extremely ugly. The following is good:                     *
 *          long k;                                                           *
 *          long N = 100;                                                     *
 *          for (k=0; k<N; ++k)                                               *
 *              do stuff                                                      *
 *      And this looks horrid:                                                *
 *          long N = 100;                                                     *
 *          for (long k=0; k<N; ++k)                                          *
 *              do stuff                                                      *
 *  5.) If a number is needed in a for loop, even once, declare or set it as  *
 *      a macro, rather than having the number verbatim in the loop. If the   *
 *      the number needs to change it's easier to keep track this way.        *
 ******************************************************************************
 *  Author:     Ryan Maguire, Wellesley College                               *
 *  Date:       June 21, 2019                                                 *
 ******************************************************************************
 *                                History                                     *
 ******************************************************************************
 *  2019/06/21 (Ryan Maguire):                                                *
 *      Initial commit.                                                       *
 *  2020/07/28 (Ryan Maguire):                                                *
 *      Clarified comments, fixed bug in error checks.                        *
 *  2020/08/22 (Ryan Maguire):                                                *
 *      Added FFT routine.                                                    *
 *  2020/09/06 (Ryan Maguire):                                                *
 *      Removed FFTW dependence. Replaced with new rss_ringoccs FFT routine.  *
 ******************************************************************************/
#include <stdlib.h>
#include <libtmpl/include/tmpl.h>
#include <stdio.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      DiffractionCorrectionFresnel                                          *
 *  Purpose:                                                                  *
 *      Compute the Fresnel transform using the classic Fresnel quadratic     *
 *      approximation to the Fresnel kernel.                                  *
 *  Arguments:                                                                *
 *      dlp (DLPObj *):                                                       *
 *          An instance of the DLPObj structure defined in                    *
 *          rss_ringoccs_diffraction_correction.h. This contains all of the   *
 *          necessary data for diffraction correction, including the geometry *
 *          of the occultation and actual power and phase data.               *
 *  Output:                                                                   *
 *      Nothing (void):                                                       *
 *          This is a void function, so no actual output is provided. However *
 *          the T_out pointer within the dlp structure will be changed at the *
 *          end, containing the diffraction corrected data.                   *
 *  Notes:                                                                    *
 *      1.) This code uses the Fresnel approximation which has been known to  *
 *          fail for several different occultations, especially ones of very  *
 *          low angle (small B values). Take this into consideration when     *
 *          performing any analysis.                                          *
 *      2.) While this may be inaccurate for certain occultations, it is      *
 *          immensely fast, capable of processing the entire Rev007 E         *
 *          occultation accurately in less than a second at 1km resolution.   *
 ******************************************************************************/
void rssringoccs_Diffraction_Correction_Fresnel(rssringoccs_TAUObj *tau)
{
    /*  m and n used for indexing, nw_pts is number of points in window.      */
    size_t m, n, nw_pts, center;

    /*  w_init is window width (km), dx and two_dx are sample spacing (km).   */
    double w_init, dx, two_dx;

    /*  Pointers for the independent variable and the window function.        */
    double *x_arr;
    double *w_func;
    double fwd_factor;

    /*  Declare the window function pointer and allocate memory for it. The   *
     *  type rssringoccs_window_func was declared at the start of this file.  *
     *  Be sure to free this at the end!                                      */
    rssringoccs_Window_Function fw = tau->window_func;

    void (*FresT)(
        rssringoccs_TAUObj *, const double *, const double *, size_t, size_t
    );

    /*  This should remain at false.                                          */
    tau->error_occurred = tmpl_False;

    if (tau->use_norm)
        FresT = rssringoccs_Fresnel_Transform_Norm;
    else
        FresT = rssringoccs_Fresnel_Transform;

    if (tau->use_fwd)
        fwd_factor = -1.0;
    else
        fwd_factor = 1.0;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Data(tau);

    if (tau->error_occurred)
        return;

    /*  Retrieve the starting point from the TAUObj instance.                 */
    center = tau->start;

    /*  Compute some extra necessary variables. At this point it is assumed   *
     *  that tau->w_km_vals, tau->rho_km_vals, and others are pointers, most  *
     *  likely created with malloc or calloc, that point to a memory block    *
     *  that is tau->arr_size in size, that tau->start >= 0, and that         *
     *  tau->start+tau->n_used <= tau->arr_size. No error checks for this     *
     *  are performed here, but rather the caller of this function has that   *
     *  responsibility. Such checks are performed in the                      *
     *  DiffractionCorrection Python class, so if you're only using that then *
     *  there's no problem. If not, this next step may cause a segmentation   *
     *  fault.                                                                */
    w_init = tau->w_km_vals[center];

    /*  It is also assumed these pointers have at least two elements of       *
     *  doubles being pointed to. Again, DiffractionCorrection checks this.   *
     *  Hence dlp->rho_km_vals[center] and dlp->rho_km_vals[center+1] should  *
     *  be initialized. If not, you will get a segmentation fault.            */
    dx = tau->dx_km;

    /*  And now, some more variables we'll be using.                          */
    two_dx = 2.0*dx;
    nw_pts = (size_t)(w_init / two_dx) + 1UL;

    /* Check to ensure you have enough data to the left.                      */
    rssringoccs_Tau_Check_Data_Range(tau);

    if (tau->error_occurred)
        return;

    /*  Reserve some memory for two arrays, the ring radius and the window    *
     *  function. This will need to be reallocated later if the window width  *
     *  changes by more than two_dx. Per C99 recommendations, we do not cast  *
     *  malloc since void will safely be promoted to a double pointer.        */
    x_arr  = malloc(sizeof(*x_arr)  * nw_pts);
    w_func = malloc(sizeof(*w_func) * nw_pts);

    /*  Check that malloc was successfull then pass the x_arr array           *
     *  (ring radius) to the void function reset_window. This alters x_arr so *
     *  that it's values range from -W/2 to zero, W begin the window width.   */
    if (!(x_arr))
    {
        tau->error_occurred = tmpl_True;
        tau->error_message = tmpl_strdup(
            "\n\rError Encountered: rss_ringoccs\n\n"
            "\r\tDiffractionCorrectionFresnel\n\n"
            "\rMalloc failed and returned NULL for x_arr. Returning.\n\n"
        );
        return;
    }

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

    rssringoccs_Tau_Reset_Window(x_arr, w_func, dx, w_init, nw_pts, fw);

    /* Compute Window Functions, and compute pi/2 * x^2                       */
    for (m = 0; m < nw_pts; ++m)
    {
        /*  The independent variable is pi/2 * ((rho-rho0)/F)^2. Compute      *
         *  part of this. The 1/F^2 part is introduced later.                 */
        x_arr[m] *= tmpl_Pi_By_Two*x_arr[m];

        /*  Use the fwd_factor to computer forward or inverse transform.      */
        x_arr[m] *= fwd_factor;
    }

    /*  Compute the Fresnel transform across the input data.                  */
    for (m=0; m<=tau->n_used; ++m)
    {
        /*  If the window width has deviated more the 2*dx, reset values.     */
        if (tmpl_Double_Abs(w_init - tau->w_km_vals[center]) >= two_dx)
        {
            /* Reset w_init and recompute window function.                    */
            w_init = tau->w_km_vals[center];
            nw_pts = ((size_t)(w_init / two_dx)) + 1UL;

            /*  Reallocate memory, since the sizes of the arrays changed.     */
            w_func = realloc(w_func, sizeof(double)*nw_pts);
            x_arr  = realloc(x_arr, sizeof(double)*nw_pts);

            /*  Reset the x_arr array to range between -W/2 and zero.         */
            rssringoccs_Tau_Reset_Window(x_arr, w_func, dx, w_init, nw_pts, fw);

            /* Compute Window Functions, and compute pi/2 * x^2.              */
            for(n=0; n<nw_pts; ++n)
            {
                x_arr[n] *= tmpl_Pi_By_Two*x_arr[n];

                /*  Again, if forward calculation is set, negate x_arr.       */
                x_arr[n] *= fwd_factor;
            }
        }

        /*  Compute the Fresnel Transform about the current point.            */
        FresT(tau, x_arr, w_func, nw_pts, center);

        /*  Move the pointers to the next point.                              */
        center += 1;
    }

    /*  Free the variables allocated by malloc.                               */
    free(x_arr);
    free(w_func);
}
