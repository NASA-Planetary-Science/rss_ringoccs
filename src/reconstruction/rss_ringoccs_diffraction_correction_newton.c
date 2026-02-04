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
 ******************************************************************************/
#include <libtmpl/include/compat/tmpl_cast.h>
#include <libtmpl/include/compat/tmpl_free.h>
#include <libtmpl/include/compat/tmpl_malloc.h>
#include <libtmpl/include/tmpl_bool.h>
#include <libtmpl/include/tmpl_math.h>
#include <rss_ringoccs/include/rss_ringoccs_tau.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>
#include <stddef.h>

/*  Performs the Newton transform on a Tau object.                            */
void rssringoccs_Diffraction_Correction_Newton(rssringoccs_TAUObj * const tau)
{
    /*  Variable for indexing the inner most loop.                            */
    size_t n;

    /*  Check that tau isn't NULL before trying to access its members.        */
    if (!tau)
        return;

    /*  Check that the pointers to the data are not NULL.                     */
    rssringoccs_Tau_Check_Core_Data(tau);

    /* Check to ensure you have enough data to process the transform.         */
    rssringoccs_Tau_Check_Data_Range(tau);

    /*  The previous functions set the error_occurred Boolean to true on      *
     *  error. Do not proceed if an error occurred, check for this.           */
    if (tau->error_occurred)
        return;

    /*  Scope for declaring const variables. C89 does not permitting mixing   *
     *  code with declarations, and we needed to run error checks above, so   *
     *  we are unable to declare these variables with the const qualifier     *
     *  without creating a new scope (made by using curly braces, {}).        *
     *  Since we want to enable parallelization, telling the compiler that    *
     *  these variables are constant and different threads may safely read    *
     *  from them is desirable. For portability with compilers that do not    *
     *  support recent standards, start a new scope for these declarations.   */
    {
        /*  The rho_km_vals array is required to have equally spaced samples. *
         *  The displacement between bins can be computed from the difference *
         *  of any two successive points.                                     */
        const double dx = tau->rho_km_vals[1] - tau->rho_km_vals[0];
        const double two_dx = 2.0 * dx;

        /*  When the window width is large, or the sample spacing is very     *
         *  fine, repeated calls to realloc, which resizes the window array,  *
         *  start to bottleneck the computation. It is faster to compute the  *
         *  maximum width required for the entire computation and then make   *
         *  one call to malloc. These variables are for the largest window    *
         *  width required.                                                   */
        const double * const w_array = tau->w_km_vals + tau->start;
        const double w_max = tmpl_Double_Array_Max(w_array, tau->n_used);

        /*  We need an odd number of points for the window. N to the right of *
         *  the center, N to the left, and the center itself. N is obtained   *
         *  by casting w_max / two_dx to an integer.                          */
        const size_t nw_max = (TMPL_CAST(w_max / two_dx, size_t) << 1) + 1;

        /*  Select the correct function pointer.                              */
        const rssringoccs_FresnelNewtonTransform newton_transform =
            rssringoccs_Tau_Select_Newton_Transform(tau);

        /*  The selector function checks if the psinum in the Tau object does *
         *  not represent a Newton transform. It will set the error_occurred  *
         *  Boolean should this take place. Check.                            */
        if (tau->error_occurred)
            return;

        /*  At the time of this writing, not every Newton transform has been  *
         *  implemented (there are many!). If the selector function returned  *
         *  NULL, the user requested a psitype that is not available.         */
        if (!newton_transform)
        {
            tau->error_occurred = tmpl_True;
            tau->error_message =
                "\n\rError Encountered: rss_ringoccs\n"
                "\r\trssringoccs_Diffraction_Correction_Newton\n\n"
                "\rrssringoccs_Tau_Select_Newton_Transform returned NULL.\n"
                "\rRequested psitype not implemented yet.\n";

            return;
        }

        /*  Begin a new scope for variables that each thread should get its   *
         *  own copy of. If OpenMP is not enabled, this has no effect, but if *
         *  it is it allows us to parallelize the main for-loop safely.       */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            /*  w_init is the initial window width. We do not need to         *
             *  repeatedly update the window function, we only need to        *
             *  recompute once the width grows by at least two bins (one on   *
             *  the left and one on the right). That is, once                 *
             *  |w_km_vals[center] - w_init| > 2 * dx is true.                */
            double w_init = tau->w_km_vals[tau->start];

            /*  The windows must contain an odd number of samples, half to    *
             *  the left of the center, half to the right, and the center.    */
            size_t nw_pts = (TMPL_CAST(w_init / two_dx, size_t) << 1) + 1;

            /*  Pointer to an array for the tapering function. w_func[n]      *
             *  corresponds to the evaluation of w(x[n]) at the nth point in  *
             *  the window. w is often some combination of Bessel functions,  *
             *  it is beneficial to precompute the tapering function once     *
             *  and use that.                                                 */
            double *w_func = TMPL_MALLOC(double, nw_max);

            /*  Malloc returns NULL on failure. Check for this.               */
            if (!w_func)
            {
                tau->error_occurred = tmpl_True;
                tau->error_message =
                    "\n\rError Encountered: rss_ringoccs\n"
                    "\r\trssringoccs_Diffraction_Correction_Newton\n\n"
                    "\rmalloc returned NULL for w_func.\n";

                /*  In-case OpenMP support is enabled, we can not break out   *
                 *  of this scope and simply return to the caller. Jump ahead *
                 *  to the cleanup stage, which is after the for-loop.        */
                goto CLEANUP;
            }

            /*  Other threads should abort if one of the others failed.       */
            if (tau->error_occurred)
                goto CLEANUP;

            /*  Compute the values in the window function. Note, w_func       *
             *  almost surely has more than nw_pts allocated to it, it has    *
             *  nw_max allocated. This function sets the first nw_pts in the  *
             *  array, the remaining points are uninitialized. This makes no  *
             *  difference, all of the routines used take in the nw_pts       *
             *  variable as a parameter and do not read beyond that length.   */
            rssringoccs_Tau_Compute_Window(tau, w_func, nw_pts, tau->start);

            /*  Loop through all of the points in the processing range and    *
             *  run the transform (either forward or inverse) on the data.    *
             *  We have initialized the w_func array, and the corresponding   *
             *  variables w_init and nw_pts, inside a scope that starts with  *
             *  #pragma omp parallel. Each thread get's its own copy of those *
             *  variables. The functions in the inner-most for-loop are both  *
             *  thread-safe and reentrant, we may parallelize the computation.*/
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (n = 0; n < tau->n_used; ++n)
            {
                /*  The window width varies with the Fresnel scale, which     *
                 *  depends on the geometry and hence varies across the data  *
                 *  set. Compute how much the width has changed so far.       */
                const size_t center = tau->start + n;
                const double w_diff = w_init - tau->w_km_vals[center];
                const double w_err = tmpl_Double_Abs(w_diff);

                /*  If the window width changed too much, recompute w_func.   */
                if (w_err >= two_dx)
                {
                    /* Reset w_init and recompute window function.            */
                    w_init = tau->w_km_vals[center];
                    nw_pts = (TMPL_CAST(w_init / two_dx, size_t) << 1) + 1;

                    /*  Compute the window function for the new width.        */
                    rssringoccs_Tau_Compute_Window(tau, w_func, nw_pts, center);
                }

                /*  Compute the Newton transform about the current point.     */
                newton_transform(tau, w_func, nw_pts, center);
            }
            /*  End of for-loop for n.                                        */

            /*  Free variables allocated by malloc.                           */
            CLEANUP:
                TMPL_FREE(w_func);
        }
        /*  End of #pragma omp parallel scope.                                */
    }
    /*  End of main parameter declaration scope.                              */
}
/*  End of rssringoccs_Diffraction_Correction_Newton.                         */