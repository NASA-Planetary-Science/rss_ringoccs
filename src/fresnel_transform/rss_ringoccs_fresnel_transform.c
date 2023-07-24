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
#include <libtmpl/include/tmpl.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_transform.h>

/******************************************************************************
 *  Function:                                                                 *
 *      Fresnel_Transform_Double / Fresnel_Transform_Norm_Double              *
 *  Purpose:                                                                  *
 *      Perform diffraction correction on diffraction limited data using the  *
 *      classic Fresnel quadratic approximation to the Fresnel kernel.        *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as pi/2 (rho-rho0)^2, where rho is the ring radius of     *
 *          the point being reconstructed, and rho0 is the dummy variable of  *
 *          integration which varies from rho-W/2 to zero, W being the        *
 *          window width. This should point to n_pts number of elements.      *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. This error check check is perform  *
 *          when the function is called in _diffraction_functions.c, as well  *
 *          as in the Python DiffractionCorrection class found in             *
 *          diffraction_correction.py.                                        *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr. This should *
 *          point to n_pts number of elements.                                *
 *      F (double):                                                           *
 *          The Fresnel scale at rho.                                         *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      n_pts (long):                                                         *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *          with only half of the window. This saves a lot of computation.    *
 *      T_in_steps (long):                                                    *
 *          The number of steps in memory of the nth point to the (n+1)th     *
 *          point in the T_in pointer.                                        *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile.                                *
 ******************************************************************************/
void
rssringoccs_Fresnel_Transform(rssringoccs_TAUObj *tau,
                              const double *x_arr,
                              const double *w_func,
                              size_t n_pts,
                              size_t center)
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    size_t m, n;

    /*  rcpr_F and rcpr_F2 are the reciprocal of the Fresnel scale, and the   *
     *  square of this. x is used as the argument of the Fresnel kernel.      */
    double x, rcpr_F, rcpr_F2, factor;

    /*  exp_negative_ix is used for the Fresnel kernel.                       */
    tmpl_ComplexDouble exp_negative_ix, integrand, arg;

    /*  Start with the central point in the Riemann sum. This is center of    *
     *  window function. That is, where w_func = 1. This is just T_in at      *
     *  the central point. This also initializes T_out.                       */
    tau->T_out[center] = tau->T_in[center];

    /*  From symmetry we need only compute -W/2 to zero, so start at -n_pts.  */
    n = n_pts;

    /*  Division is more expensive than multiplication, so store the          *
     *  reciprocal of F as a variable and compute with that.                  */
    rcpr_F = 1.0 / tau->F_km_vals[center];
    rcpr_F2 = rcpr_F*rcpr_F;
    factor = 0.5*tau->dx_km*rcpr_F;

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        x = x_arr[m]*rcpr_F2;

        /*  Use Euler's Theorem to compute exp(-ix). Scale by window function.*/
        exp_negative_ix = tmpl_CDouble_Polar(w_func[m], -x);

        /*  Take advantage of the symmetry of the quadratic approximation.    *
         *  This cuts the number of computations roughly in half. If the T_in *
         *  pointer does not contain at least 2*n_pts+1 points, n_pts to the  *
         *  left and n_pts to the right of the center, then this will create  *
         *  a segmentation fault, crashing the program.                       */
        integrand = tmpl_CDouble_Add(tau->T_in[center - n],
                                     tau->T_in[center + n]);
        integrand = tmpl_CDouble_Multiply(exp_negative_ix, integrand);
        tmpl_CDouble_AddTo(&tau->T_out[center], &integrand);
        n--;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    arg = tmpl_CDouble_Rect(factor, factor);
    tau->T_out[center] = tmpl_CDouble_Multiply(arg, tau->T_out[center]);
}
