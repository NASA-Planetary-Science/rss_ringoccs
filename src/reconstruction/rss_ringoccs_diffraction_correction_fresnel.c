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
 *      The inverse Fresnel transform is given by:                            *
 *                                                                            *
 *                           infinity                                         *
 *                              -                                             *
 *                     1 + i   | |                                            *
 *          T_out(r) = -----   |   T_in(r0) w(r - r0) exp(-i psi(r, r0)) dr0  *
 *                       2F  | |                                              *
 *                            -                                               *
 *                          -infinity                                         *
 *                                                                            *
 *      where T_in is the diffracted data, w is the window function, r is     *
 *      the ring intercept point, and r0 is a dummy variable of integration.  *
 *      psi is the Fresnel Kernel, and exp is the exponential function. The   *
 *      scale factor F is the Fresnel scale which is dependent on the         *
 *      geometry of the occultation.                                          *
 *                                                                            *
 *      In ideal scenarios, the Fresnel kernel may be approximated by a       *
 *      simple quadratic:                                                     *
 *                                                                            *
 *                     -      - 2                                             *
 *                pi  | r - r0 |                                              *
 *          psi = --- | ------ |                                              *
 *                 2  |    F   |                                              *
 *                     -      -                                               *
 *                                                                            *
 *      The above integral is then computed via a Riemann sum using this new  *
 *      expression for psi.                                                   *
 *                                                                            *
 *      As the resolution gets too coarse, say 10 km or larger, the window    *
 *      width quickly shrinks to zero and the integral will be approximately  *
 *      zero. To account for this, we normalize the integral by the window    *
 *      width. The normalization is defined as follows:                       *
 *                                                                            *
 *                    |  infinity                 |                           *
 *                    |     -                     |                           *
 *                    |    | |                    |                           *
 *                    |    |    exp(-i psi(x)) dx |                           *
 *                    |  | |                      |                           *
 *                    |   -                       |                           *
 *                    | -infinity                 |                           *
 *          norm =  ---------------------------------                         *
 *                  |   W / 2                       |                         *
 *                  |     -                         |                         *
 *                  |    | |                        |                         *
 *                  |    |    w(x)exp(-i psi(x)) dx |                         *
 *                  |  | |                          |                         *
 *                  |   -                           |                         *
 *                  | -W / 2                        |                         *
 *                                                                            *
 *      This has the effect of making the integral in free-space regions,     *
 *      which are regions that are not affected by diffraction, evaluate to   *
 *      one, regardless of what (positive) resolution is chosen.              *
 *  Notes:                                                                    *
 *      1.) This code uses the Fresnel approximation which has been known to  *
 *          fail for several different occultations, especially ones of very  *
 *          low angle (small B values). Take this into consideration when     *
 *          performing any analysis.                                          *
 *      2.) While this may be inaccurate for certain occultations, it is      *
 *          immensely fast, capable of processing the entire Rev007 E         *
 *          occultation accurately in less than a second at 1km resolution.   *
 *  References:                                                               *
 *      1.) Maguire, R., French, R. (2025)                                    *
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
 *  2.) compat/tmpl_cast.h:                                                   *
 *          Macros for casting with compatibility for both C and C++.         *
 *  3.) tmpl_free.h:                                                          *
 *          Provides a macro for free with C vs. C++ compatibility.           *
 *  4.) tmpl_malloc.h:                                                        *
 *          Provides a macro for malloc with C vs. C++ compatibility.         *
 *  5.) tmpl_math_constants.h:                                                *
 *          Contains common mathematical constants, such as multiples of pi.  *
 *  6.) rss_ringoccs_tau.h:                                                   *
 *          Header file where the rssringoccs_TAUObj type is provided.        *
 *  7.) rss_ringoccs_fresnel_transform.h:                                     *
 *          Header file where the Riemann sums for various types of Fresnel   *
 *          transforms are provided. These routines perform the inner most    *
 *          for-loop in the Fresnel transforms.                               *
 *  8.) rss_ringoccs_reconstruction.h:                                        *
 *          Header file with the function prototype.                          *
 *  9.) stddef.h:                                                             *
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

/*  rssringoccs_TAUObj typedef given here.                                    */
#include <rss_ringoccs/include/rss_ringoccs_tau.h>

/*  The Riemann summation routines are found here.                            */
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/*  Function prototype provided here.                                         */
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/*  size_t provided here, data type used for indexing arrays.                 */
#include <stddef.h>

/*  The inner most for-loop is a Riemann sum for the Fresnel transform. The   *
 *  user has the option to normalize this by the window width. There are      *
 *  thus two transforms: without normalization, and with normalization.       */
static const rssringoccs_FresnelTransform
rssringoccs_fresnel_transform_list[2] = {
    rssringoccs_Fresnel_Transform,
    rssringoccs_Fresnel_Transform_Norm
};

/*  Performs the Fresnel transform on the data contained in tau.              */
void rssringoccs_Diffraction_Correction_Fresnel(rssringoccs_TAUObj * const tau)
{
    /*  m and n used for indexing, n_pts is number of points in the arrays,   *
     *  and center is the index for the center of the window of integration.  */
    size_t m, n, n_pts, center;

    /*  w_init is window width (km), two_dx is double the sample spacing (km).*/
    double w_init, two_dx;

    /*  Pointers for the independent variable and the window function.        */
    double *x_arr = NULL;
    double *w_func = NULL;

    /*  The user has two options for transforms. We'll set this later.        */
    rssringoccs_FresnelTransform fresnel_transform;

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

    /*  Select the desired transform. For the Fresnel method there are only   *
     *  two options: with or without normalization. tau contains a Boolean,   *
     *  "use_norm", which can be used to index the table above. Use this.     */
    fresnel_transform = rssringoccs_fresnel_transform_list[tau->use_norm];

    /*  Retrieve the starting point from the TAUObj struct.                   */
    center = tau->start;

    /*  Compute necessary data for the start of the inversion.                */
    w_init = tau->w_km_vals[center];
    two_dx = 2.0 * tau->dx_km;
    n_pts = TMPL_CAST(w_init / two_dx, size_t) + 1;

    /*  Allocate memory for the independent variable and the window function. */
    x_arr = TMPL_MALLOC(double, n_pts);
    w_func = TMPL_MALLOC(double, n_pts);

    /*  Check if malloc failed. It returns NULL on failure.                   */
    if (!x_arr || !w_func)
    {
        /*  If we get here, malloc failed. Abort with error.                  */
        tau->error_occurred = tmpl_True;
        tau->error_message =
            "\n\rError Encountered: rss_ringoccs\n"
            "\r\trssringoccs_Diffraction_Correction_Fresnel\n\n"
            "\rMalloc failed and returned NULL. Aborting.\n\n";

        /*  It is possible that malloc succeeded for one variable and not the *
         *  other. We initialized both pointers to NULL at the start of this  *
         *  function, so we can safely free memory using the following macro. */
        TMPL_FREE(x_arr);
        TMPL_FREE(w_func);

        return;
    }

    /*  Initialize the window array and the independent variable.             */
    rssringoccs_Tau_Reset_Window(
        tau,        /*  Tau object containing the window function.            */
        x_arr,      /*  The independent variable, r[n] - r[center].           */
        w_func,     /*  The window array as a function of x_arr.              */
        n_pts,      /*  Number of points in the x_arr array.                  */
        center      /*  Index for the center of the window.                   */
    );

    /*  We have computed the window function and the independent variable x,  *
     *  which is (r - r0). We need (pi/2) (r - r0)^2. The 1 / F^2 factor      *
     *  is introduced later inside the fresnel_transform function.            */
    for (m = 0; m < n_pts; ++m)
        x_arr[m] *= tmpl_Double_Pi_By_Two * x_arr[m];

    /*  Compute the Fresnel transform across the input data.                  */
    for (n = 0; n <= tau->n_used; ++n)
    {
        /*  Check if we need to resize the window. This happens once          *
         *  |w - w0| > 2 * dx occurs, where w is the window width of the      *
         *  current index, and w0 is the width of the window being used.      */
        const int resize = rssringoccs_Tau_Resize_Half_Window(
            tau, &x_arr, &w_func, w_init, two_dx, center
        );

        /*  If we did need to resize, there are a few things that could have  *
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
                return;

            /*  Reset the threshold value for the window width to the current *
             *  window size. We will update again once the required window    *
             *  grows beyond 2 * dx the size of the current window.           */
            w_init = tau->w_km_vals[center];

            /*  Similarly, reset the number of points in the window.          */
            n_pts = TMPL_CAST(w_init / two_dx, size_t) + 1;

            /*  The rssringoccs_Tau_Resize_Half_Window function calls the     *
             *  rssringoccs_Tau_Reset_Window routine, which recomputes the    *
             *  expression (r - r0) across the new window. We need the        *
             *  expression (pi/2) (r - r0)^2. Compute this.                   */
            for (m = 0; m < n_pts; ++m)
                x_arr[m] *= tmpl_Double_Pi_By_Two * x_arr[m];
        }

        /*  Compute the Fresnel Transform about the current point.            */
        fresnel_transform(tau, x_arr, w_func, n_pts, center);

        /*  Move the pointers to the next point.                              */
        center += 1;
    }

    /*  Free the variables allocated by malloc.                               */
    TMPL_FREE(x_arr);
    TMPL_FREE(w_func);
}
/*  End of rssringoccs_Diffraction_Correction_Fresnel.                        */
