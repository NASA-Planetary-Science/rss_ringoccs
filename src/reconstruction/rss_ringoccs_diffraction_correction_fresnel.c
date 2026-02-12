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
 *                rss_ringoccs_diffraction_correction_fresnel                 *
 ******************************************************************************
 *  Purpose:                                                                  *
 *      Uses the quadratic Fresnel transform for diffraction correction.      *
 ******************************************************************************
 *                             DEFINED FUNCTIONS                              *
 ******************************************************************************
 *  Function Name:                                                            *
 *      rssringoccs_Diffraction_Correction_Fresnel                            *
 *  Purpose:                                                                  *
 *      Performs diffraction correction using the Fresnel approximation.      *
 *  Arguments:                                                                *
 *      tau (rssringoccs_TAUObj * const):                                     *
 *          The geometry and diffraction data for the reconstruction. The     *
 *          output reconstruction is stored in the T_out member of tau.       *
 *  Output:                                                                   *
 *      None (void).                                                          *
 *  Called Functions:                                                         *
 *      None.                                                                 *
 *  Method:                                                                   *
 *      The cylindrical Fresnel transform is given by:                        *
 *                                                                            *
 *                             inf   2 pi                                     *
 *                              -     -                                       *
 *          ^        sin(B)    | |   | |      exp(i  psi)                     *
 *          T(r ) = --------   |     |   T(r) ----------- r dr dphi           *
 *             0    i lambda | |   | |        | R - rho |                     *
 *                            -     -                                         *
 *                            0     0                                         *
 *                                                                            *
 *      where T_hat is the diffracted data, T is the complex transmittance of *
 *      the aperture, rho = (r, phi) is the dummy variable of integration     *
 *      across the aperture (which lies in the z = 0 plane),                  *
 *      rho0 = (r0, phi0) is the point of interest, which is the point in     *
 *      the aperture plane where the line of sight from the observer          *
 *      intersects, R = (x, y, z) is the position vector for the observer, B  *
 *      is the "opening angle," which is the angle the line of sight from the *
 *      observer makes with the z = 0 plane, psi is the Fresnel phase         *
 *      (which varies with r, r0, phi, phi0, and B), and lambda is the        *
 *      wavelength of the incident wave.                                      *
 *                                                                            *
 *      The one-dimensional transform is obtained by applying the stationary  *
 *      phase approximation to the Fresnel phase psi with respect to the      *
 *      (dummy variable of integration) azimuth angle phi. If phi_s is the    *
 *      solution to d psi / d phi = 0, if rho_s = (r, phi_s), and if psi_s is *
 *      the evaluation of psi at the tuple (r, phi_s, r0, phi0), then the     *
 *      stationary phase approximation produces:                              *
 *                                                                            *
 *            2 pi                                                            *
 *             -                          _________                           *
 *            | | exp(i  psi)            /  2 pi    exp(i (psi_s - pi/4))     *
 *            |   ----------- dphi ~=   / --------- ---------------------     *
 *          | |   | R - rho |         \/  |psi_s''|     | R - rho_s |         *
 *           -                                                                *
 *           0                                                                *
 *                                                                            *
 *      The one-dimensional transform is then:                                *
 *                                                                            *
 *                     inf                                                    *
 *                      -           _________                                 *
 *          ^          | |         /  2 pi    exp(i (psi_s - pi/4))           *
 *          T(r ) ~=   |   T(r)   / --------- --------------------- r dr      *
 *             0     | |        \/  |psi ''|      | R - rho |                 *
 *                    -                 s                  s                  *
 *                    0                                                       *
 *                                                                            *
 *      In ideal scenarios, we may apply the following approximations:        *
 *                                                                            *
 *               r              1        1                                    *
 *          ----------- ~= ----------- = -                                    *
 *          | R - rho |    | R - rho |   D                                    *
 *                   s              0                                         *
 *                                                                            *
 *              _________                                                     *
 *             /  2 pi            1                                           *
 *            / --------- ~= ---------                                        *
 *          \/  |psi_s''|      ___                                            *
 *                           \/ 2   F                                         *
 *                                                                            *
 *                      -      - 2                                            *
 *                     | r - r  |                                             *
 *                 pi  |      0 |                                             *
 *          psi ~= --- | ------ | = psi                                       *
 *             s    2  |    F   |      Q                                      *
 *                      -      -                                              *
 *                                                                            *
 *      where D is the distance from the observer (R) to the point of         *
 *      interest in the plane (rho0), psi_Q is the quadratic Fresnel phase,   *
 *      and F is the "Fresnel scale," given by:                               *
 *                                                                            *
 *                                                                            *
 *                       ______________________________                       *
 *                      /                   2         2                       *
 *                     /          1 - cos(B) sin(phi )                        *
 *                    / lambda D                    0                         *
 *          F = \    /  -------- ---------------------                        *
 *               \  /      2                  2                               *
 *                \/                    sin(B)                                *
 *                                                                            *
 *                                                                            *
 *      This combination of approximations produces the classic Fresnel       *
 *      quadratic transform:                                                  *
 *                                                                            *
 *                                         -                 -                *
 *                           inf          |        -      - 2 |               *
 *                            -           |       | r - r  |  |               *
 *          ^        1 - i   | |          |   pi  |      0 |  |               *
 *          T(r ) ~= -----   |   T(r) exp | i --- | ------ |  | dr            *
 *             0      2 F  | |            |    2  |    F   |  |               *
 *                          -             |        -      -   |               *
 *                          0              -                 -                *
 *                                                                            *
 *      In the case where B = pi / 2, or if the angle phi0 is fixed, then the *
 *      Fresnel scale becomes independent of rho0 and the resulting           *
 *      transform is a convolution:                                           *
 *                                                                            *
 *                          -                 -                               *
 *            inf          |        -      - 2 |                              *
 *             -           |       | r - r  |  |                              *
 *            | |          |   pi  |      0 |  |                              *
 *            |   T(r) exp | i --- | ------ |  | dr = T * exp(i psi )         *
 *          | |            |    2  |    F   |  |                   Q          *
 *           -             |        -      -   |                              *
 *           0              -                 -                               *
 *                                                                            *
 *      where * denotes the convolution. Using the convolution theorem, we    *
 *      apply the Fourier transform to this allowing us to solve for T in     *
 *      terms of T_hat. That is, we may solve this inverse problem exactly:   *
 *                                                                            *
 *                                         -                  -               *
 *                          inf           |         -      - 2 |              *
 *                           -            |        | r - r  |  |              *
 *                  1 + i   | | ^         |    pi  |      0 |  |              *
 *          T(r) ~= -----   |   T(r ) exp | -i --- | ------ |  | dr           *
 *                   2 F  | |      0      |     2  |    F   |  |   0          *
 *                         -              |         -      -   |              *
 *                         0               -                  -               *
 *                                                                            *
 *      This function uses this inverse formula to reconstruct T from T_hat.  *
 *      Since the window of integration is finite, a tapering function is     *
 *      introduced to removed unwanted artifacts. The transform becomes:      *
 *                                                                            *
 *                                                  -                  -      *
 *                       r + W/2                   |         -      - 2 |     *
 *                          -                      |        | r - r  |  |     *
 *                 1 + i   | | ^                   |    pi  |      0 |  |     *
 *          T(r) = -----   |   T(r ) w(r - r ) exp | -i --- | ------ |  | dr  *
 *                  2 F  | |      0         0      |     2  |    F   |  |   0 *
 *                        -                        |         -      -   |     *
 *                     r - W/2                      -                  -      *
 *                                                                            *
 *      As the resolution gets too coarse, say 10 km or larger, the window    *
 *      width quickly shrinks to zero and the integral will be approximately  *
 *      zero. To account for this, the integral may be optionally normalized  *
 *      by the window width. The is defined via:                              *
 *                                                                            *
 *                    |    inf                    |                           *
 *                    |     -                     |                           *
 *                    |    | |                    |                           *
 *                    |    |    exp(-i psi(x)) dx |                           *
 *                    |  | |                      |                           *
 *                    |   -                       |                           *
 *                    | -inf                      |                           *
 *          norm =  ---------------------------------                         *
 *                  |   W / 2                       |                         *
 *                  |     -                         |                         *
 *                  |    | |                        |                         *
 *                  |    |   w(x) exp(-i psi(x)) dx |                         *
 *                  |  | |                          |                         *
 *                  |   -                           |                         *
 *                  | -W / 2                        |                         *
 *                                                                            *
 *      This has the effect of making the integral in free-space regions,     *
 *      which are regions that are not affected by diffraction, evaluate to   *
 *      one, regardless of what (positive) resolution is chosen. Set the      *
 *      normalize Boolean to true to enable this.                             *
 *  Notes:                                                                    *
 *      1.) This code uses the Fresnel approximation which has been known to  *
 *          fail for several different occultations, especially ones of very  *
 *          low angle (small B values). Take this into consideration when     *
 *          performing any analysis.                                          *
 *                                                                            *
 *      2.) While this may be inaccurate for certain occultations, it is      *
 *          immensely fast, capable of processing the entire Rev007 E         *
 *          occultation accurately in a fraction of a second with 250 meter   *
 *          sampling and 1000 meter resolution.                               *
 *                                                                            *
 *      3.) The inner for-loop may be optionally parallelized. To enable      *
 *          this, ensure that rss_ringoccs was compiled with OpenMP support.  *
 *  References:                                                               *
 *      1.) Maguire, R., French, R. (2026)                                    *
 *          Applications of Legendre Polynomials for Fresnel Inversion        *
 *              and Occultation Observations                                  *
 *                                                                            *
 *          A full derivation of the Fresnel approximation is given here.     *
 *                                                                            *
 *      2.) Marouf, E., Tyler, G., Rosen, P. (June 1986)                      *
 *          Profiling Saturn's Rings by Radio Occultation                     *
 *          Icarus Vol. 68, Pages 120-166.                                    *
 *                                                                            *
 *          This paper describes the theory of diffraction as applied to      *
 *          planetary ring systems. rss_ringoccs implements many of the       *
 *          ideas found in this article.                                      *
 *                                                                            *
 *      3.) Goodman, J. (2005)                                                *
 *          Introduction to Fourier Optics                                    *
 *          McGraw-Hill Series in Electrical and Computer Engineering.        *
 *                                                                            *
 *          Covers most of the theory behind diffraction and the application  *
 *          of Fourier analysis to optics. The Fresnel transform is given an  *
 *          in-depth treatise in this book.                                   *
 ******************************************************************************
 *                                DEPENDENCIES                                *
 ******************************************************************************
 *  1.) tmpl_bool.h:                                                          *
 *          Header file providing Booleans (True and False).                  *
 *  2.) tmpl_cast.h:                                                          *
 *          Macros for casting with compatibility for both C and C++.         *
 *  3.) tmpl_free.h:                                                          *
 *          Provides a macro for free with C vs. C++ compatibility.           *
 *  4.) tmpl_malloc.h:                                                        *
 *          Provides a macro for malloc with C vs. C++ compatibility.         *
 *  5.) tmpl_math_constants.h:                                                *
 *          Contains common mathematical constants, such as multiples of pi.  *
 *  5.) tmpl_math.h:                                                          *
 *          Provides common mathematical routines.                            *
 *  7.) rss_ringoccs_tau.h:                                                   *
 *          Header file where the rssringoccs_TAUObj type is provided.        *
 *  8.) rss_ringoccs_fresnel_transform.h:                                     *
 *          Header file where the Riemann sums for various types of Fresnel   *
 *          transforms are provided. These routines perform the inner most    *
 *          for-loop in the Fresnel transforms.                               *
 *  9.) rss_ringoccs_reconstruction.h:                                        *
 *          Header file with the function prototype.                          *
 *  10.) stddef.h:                                                            *
 *          Standard header file providing the size_t typedef.                *
 ******************************************************************************
 *  Author:     Ryan Maguire                                                  *
 *  Date:       June 21, 2019                                                 *
 ******************************************************************************
 *                              Revision History                              *
 ******************************************************************************
 *  2019/06/21: Ryan Maguire                                                  *
 *      Initial commit.                                                       *
 *  2020/07/28: Ryan Maguire                                                  *
 *      Clarified comments, fixed bug in error checks.                        *
 *  2020/08/22: Ryan Maguire                                                  *
 *      Added FFT routine.                                                    *
 *  2020/09/06: Ryan Maguire                                                  *
 *      Removed FFTW dependence. Replaced with new rss_ringoccs FFT routine.  *
 *  2024/12/23: Ryan Maguire                                                  *
 *      Restructured the code. Added references, cleaned up comments and      *
 *      includes. FFT method moved to its own file.                           *
 *  2025/04/15: Ryan Maguire                                                  *
 *      General cleanup, changing error_message to a const char pointer, no   *
 *      longer requires a call to tmpl_String_Duplicate. Fixed some typos.    *
 *  2026/02/12: Ryan Maguire                                                  *
 *      Added OpenMP support, added more documentation.                       *
 ******************************************************************************/

/*  Booleans (True / False) provided here.                                    */
#include <libtmpl/include/tmpl_bool.h>

/*  Macros for C vs. C++ compatibility with casting data types.               */
#include <libtmpl/include/compat/tmpl_cast.h>

/*  Macros for C vs. C++ compatibility with malloc and free.                  */
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <libtmpl/include/compat/tmpl_free.h>

/*  Provides common multiples of pi, including pi / 2.                        */
#include <libtmpl/include/constants/tmpl_math_constants.h>

/*  Absolute value function, and more, provided here.                         */
#include <libtmpl/include/tmpl_math.h>

/*  rssringoccs_TAUObj typedef given here.                                    */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  The Riemann summation routines are found here.                            */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Function prototype provided here.                                         */
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  size_t provided here, data type used for indexing arrays.                 */
#include <stddef.h>

/*  Performs the Fresnel transform on the data contained in tau.              */
void rssringoccs_Diffraction_Correction_Fresnel(rssringoccs_TAUObj * const tau)
{
    /*  Variable for indexing the inner for-loops.                            */
    size_t n;

    /*  Check that the input is not NULL before attempting to access it.      */
    if (!tau)
        return;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Core_Data(tau);

    /*  Check to ensure you have enough data to process.                      */
    rssringoccs_Tau_Check_Data_Range(tau);

    /*  The previous functions set the error_occurred Boolean on error. Check.*/
    if (tau->error_occurred)
        return;

    /*  Begin a new scope for variables that each thread should get its own   *
     *  copy of. If OpenMP is not enabled, this has no effect, but if it is   *
     *  it allows us to parallelize the main for-loop safely.                 */
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        /*  Twice the sample spacing (distance between points), in km.        */
        const double two_dx = 2.0 * tau->dx_km;

        /*  When the window width is large, or the sample spacing is very     *
         *  fine, repeated calls to realloc, which resizes the window array,  *
         *  start to bottleneck the computation. It is faster to compute the  *
         *  maximum width required for the entire computation and then make   *
         *  one call to malloc. These variables are for the largest window    *
         *  width required.                                                   */
        const double * const w_array = tau->w_km_vals + tau->start;
        const double w_max = tmpl_Double_Array_Max(w_array, tau->n_used);

        /*  We set the initial window width to a negative number that is      *
         *  guaranteed to have | w_init - w[n] | > 2 * dx. This is to save on *
         *  redundant computations should the inner for-loop be parallelized. *
         *  Each thread will start with | w_init - w[n] | > 2 * dx being true *
         *  and will compute the initial window at the starting index for the *
         *  respective thread. Having each thread initialized to the starting *
         *  window width (tau->w_km_vals[tau->start]) is wasteful since some  *
         *  threads will start with larger indices and will instantly need to *
         *  recompute the window function. We set the nw_pts variable to zero *
         *  as well. Each thread will get its own copy and this will be       *
         *  initialized at the start of the inner for-loop.                   */
        double w_init = -4.0 * w_max;
        size_t nw_pts = 0;

        /*  Half the numbers of bins (or points) in the largest window. The   *
         *  windows are symmetric, with 2 * nw_pts + 1 points total. Since    *
         *  the Fresnel quadratic approximation is also symmetric (since      *
         *  x^2 = (-x)^2), we only need to compute across the left half.      */
        const size_t nw_max = TMPL_CAST(w_max / two_dx, size_t);

        /*  Allocate memory for the independent variable and window function. */
        double *x_arr = TMPL_MALLOC(double, nw_max);
        double *w_func = TMPL_MALLOC(double, nw_max);

        /*  Check if malloc failed. It returns NULL on failure. Different     *
         *  threads share access to the Tau object (if OpenMP support is      *
         *  enabled), ensure that if one thread fails, then they all do.      */
#ifdef _OPENMP
#pragma omp critical
#endif
        if (!x_arr || !w_func)
        {
            /*  If we get here, malloc failed. Abort with error.              */
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Diffraction_Correction_Fresnel\n\n"
                "\rMalloc failed and returned NULL.\n\n";
        }

        /*  In case OpenMP support is enabled, we can not break out of this   *
         *  scope and simply return to the caller. Jump ahead to the cleanup  *
         *  stage, which is after the for-loop.                               */
        if (tau->error_occurred)
            goto CLEANUP;

        /*  Compute the Fresnel transform across the input data.              */
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (n = 0; n <= tau->n_used; ++n)
        {
            /*  The window width varies with the Fresnel scale, which depends *
             *  on the geometry and hence varies across the data set. Compute *
             *  how much the width has changed so far.                        */
            const size_t center = tau->start + n;
            const double w_diff = w_init - tau->w_km_vals[center];
            const double w_err = tmpl_Double_Abs(w_diff);

            /*  If the window width changed too much, recompute w_func.       */
            if (w_err >= two_dx)
            {
                /*  Reset the threshold value for the window width to the     *
                 *  current window size. We will update again once the        *
                 *  required window grows beyond 2 * dx the size of the       *
                 *  current window.                                           */
                w_init = tau->w_km_vals[center];

                /*  Similarly, reset the number of points in the window.      */
                nw_pts = TMPL_CAST(w_init / two_dx, size_t);

                /*  Evaluate the window function across the current window    *
                 *  and reset the x array.                                    */
                rssringoccs_Tau_Reset_Fresnel_Window(
                    tau,    /*  Tau object containing the window function.    */
                    x_arr,  /*  The independent variable, r[n] - r[center].   */
                    w_func, /*  The window array as a function of x_arr.      */
                    nw_pts, /*  Number of points in the x_arr array.          */
                    center  /*  Index for the center of the window.           */
                );
            }

            /*  Compute the Fresnel Transform about the current point.        */
            rssringoccs_Fresnel_Transform(tau, x_arr, w_func, center, nw_pts);
        }

        /*  Free the variables allocated by malloc.                           */
        CLEANUP:

            /*  For each pointer, either malloc succeeded and the the pointer *
             *  points to valid data, or malloc failed and the pointer is now *
             *  NULL. In either case we may safely use the TMPL_FREE macro.   */
            TMPL_FREE(x_arr);
            TMPL_FREE(w_func);
    }
    /*  End of #pragma omp parallel scope.                                    */
}
/*  End of rssringoccs_Diffraction_Correction_Fresnel.                        */
