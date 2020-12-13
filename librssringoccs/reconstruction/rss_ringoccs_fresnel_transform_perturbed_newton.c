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
 ******************************************************************************/

#include <rss_ringoccs/include/rss_ringoccs_math.h>
#include <rss_ringoccs/include/rss_ringoccs_complex.h>
#include <rss_ringoccs/include/rss_ringoccs_fresnel_kernel.h>
#include <rss_ringoccs/include/rss_ringoccs_reconstruction.h>

/******************************************************************************
 *  Function:                                                                 *
 *      Fresnel_Transform_Newton_Double / Fresnel_Transform_Newton_Norm_Double*
 *  Purpose:                                                                  *
 *      Performs diffraction correction by using the Newton-Raphson method of *
 *      root-finding to compute the stationary value of the Fresnel kernel.   *
 *  Arguments:                                                                *
 *      x_arr (double *):                                                     *
 *          Defined as rho-rho0, where rho0 is the ring radius of the point   *
 *          being reconstructed and rho is the dummy variable of integration  *
 *          which varies from rho0-W/2 to rho+W/2, W being the window width.  *
 *          This should contain n_pts (see below) number of elements.         *
 *      phi_arr (double *):                                                   *
 *          The ring azimuth angle corresponding to rho, in radians.          *
 *          This should contain n_pts (see below) number of elements.         *
 *      T_in (char *):                                                        *
 *          The diffracted data. Must contain at least n_pts (see below)      *
 *          points to the left and right of the starting point, or a          *
 *          segmentation fault will occur. This error check check is perform  *
 *          when the function is called in _diffraction_functions.c, as well  *
 *          as in the Python DiffractionCorrection class found in             *
 *          diffraction_correction.py.                                        *
 *      w_func (double *):                                                    *
 *          The window/tapering function, as a function of x_arr.             *
 *          This should contain n_pts (see below) number of elements.         *
 *      kd (double):                                                          *
 *          The wavenumber k scaled by D (See above).                         *
 *      r (double):                                                           *
 *          The radius of the current point be reconstructed.                 *
 *      B (double):                                                           *
 *          The ring opening angle, in radians.                               *
 *      D (double):                                                           *
 *          The distance from the spacecraft to the ring-intercept point.     *
 *      EPS (double):                                                         *
 *          The allowable error in the computation of the stationary value of *
 *          the Fresnel kernel.                                               *
 *      toler (long):                                                         *
 *          The maximum number of iterations the Newton-Raphson method is     *
 *          allowed to undergo.                                               *
 *      dx (double):                                                          *
 *          The sample spacing, equivalent to x_arr[1] - x_arr[0].            *
 *      F (double):                                                           *
 *          The Fresnel scale, in the same units as D and dx.                 *
 *      n_pts (long):                                                         *
 *          Half the number of points in the window width. The symmetry of    *
 *          the quadratic approximation allows one to perform the inversion   *
 *      T_in_steps (long):                                                    *
 *          The number of steps in memory of the nth point to the (n+1)th     *
 *          point in the T_in pointer.                                        *
 *  Outputs:                                                                  *
 *      T_out (complex double):                                               *
 *          The diffraction corrected profile.                                *
 ******************************************************************************/
rssringoccs_ComplexDouble
Fresnel_Transform_Perturbed_Newton_Double(double *x_arr, double *phi_arr,
                                          rssringoccs_ComplexDouble *T_in,
                                          double *w_func, double kD,
                                          double r, double B, double D,
                                          double EPS, unsigned long toler,
                                          double dx, double F,
                                          unsigned long n_pts,
                                          unsigned long center,
                                          double perturb[5])
{
    /*  Declare all necessary variables. i and j are used for indexing.       */
    unsigned long m, n;

    /*  The Fresnel kernel and ring azimuth angle.                            */
    double psi, phi, x, poly, cos_psi, sin_psi, rcpr_F;
    rssringoccs_ComplexDouble T_out, exp_psi, integrand;

    /*  Initialize T_out and norm to zero so we can loop over later.          */
    T_out = rssringoccs_CDouble_Zero;
    rcpr_F  = 1.0/F;

    /*  Symmetry is lost without the Legendre polynomials, or Fresnel         *
     *  quadratic. Must compute everything from -W/2 to W/2.                  */
    n = center-(long)((n_pts-1)/2);

    /*  Use a Riemann Sum to approximate the Fresnel Inverse Integral.        */
    for (m = 0; m<n_pts; ++m)
    {
        /*  Factor for the polynomial perturbation.                           */
        x = (r-x_arr[m])/D;

        /*  Calculate the stationary value of psi with respect to phi.        */
        phi = Newton_Raphson_Fresnel_Psi(kD, r, x_arr[m], phi_arr[m],
                                         phi_arr[m], B, D, EPS, toler);

        /*  Compute psi and perturb by the requested polynomial.              */
        psi = rssringoccs_Double_Fresnel_Psi(kD, r, x_arr[m], phi,
                                             phi_arr[m], B, D);

        /*  Use Horner's method to compute the polynomial.                    */
        poly  = x*perturb[4]+perturb[3];
        poly  = poly*x + perturb[2];
        poly  = poly*x + perturb[1];
        poly  = poly*x + perturb[0];
        poly *= kD;
        psi += poly;

        /*  Compute the left side of exp(-ipsi) using Euler's Formula.        */
        cos_psi = rssringoccs_Double_Cos(psi);
        sin_psi = rssringoccs_Double_Sin(psi);
        exp_psi = rssringoccs_CDouble_Rect(cos_psi, -sin_psi);
        exp_psi = rssringoccs_CDouble_Multiply_Real(w_func[m], exp_psi);

        /*  Compute the transform with a Riemann sum. If the T_in pointer     *
         *  does not contain at least 2*n_pts+1 points, n_pts to the left and *
         *  right of the center, then this will create a segmentation fault.  */
        integrand = rssringoccs_CDouble_Multiply(exp_psi, T_in[m]);
        T_out     = rssringoccs_CDouble_Add(T_out, integrand);
        n += 1;
    }

    /*  Multiply result by the coefficient found in the Fresnel inverse.      */
    integrand = rssringoccs_CDouble_Rect(0.5*dx*rcpr_F, 0.5*dx*rcpr_F);
    T_out     = rssringoccs_CDouble_Multiply(integrand, T_out);
    return T_out;
}
