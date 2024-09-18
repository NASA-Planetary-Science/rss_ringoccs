"""
################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or              #
#   modify it under the terms of the GNU General Public License as published   #
#   by the Free Software Foundation, either version 3 of the License, or       #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful             #
#   but WITHOUT ANY WARRANTY# without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Purpose:                                                                   #
#       Computes the forward model of a reconstructed data set.                #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
# PDS names are all caps. Pylint doesn't like this.
# pylint: disable = invalid-name
import numpy
from . import normalize_window
from . import windows

def fresnel_forward(rho_km_vals, F_km_vals, phi_rad_vals, B_rad_vals,
                    D_km_vals, T_vals, lambda_km_vals, w_km_vals,
                    dx_km, wtype, start, n_used, norm = True,
                    fft = False, verbose = True, psitype = 'full'):
    """
        Procedure: fresnel_forward
        Purpose:   Computes the forward model of diffraction from a set of
            reconstructed data.
        Variables:
            rho_vals:        Ring radius, in kilometers.
            F_vals:          Fresnel scale, in kilometers.
            phi_rad_vals:    Ring azimuth angle, in radians.
            B_rad_vals:      Ring opening angle, in radians.
            lambda_sky_vals: Wavelength of recieved signal, in kilometers.
            D_vals:          Spacecraft-RIP distance, in kilometers.
            dx:              Sampling spacing, in kilometers.
            T_vals:          Reconstructed complex transmittance.
            w_vals:          Window width, in kilometers.
            wtype:           Window used in reconstruction, string.
            start:           Starting point of reconstructed data.
            n_used:          Number of reconstructed points.
        Keywords:
            Normalize: Parameter for normalizing the complex transmittance by
                the window function that is used. Default is True. Set to False
                to skip this feature.
        Output:
            phase_fwd_vals  - Phase of the forward model, in radians.
            T_hat_fwd_vals  - Complex transmittance of forward model.
            p_norm_fwd_vals - Normalized power of forward model, unitless.
        History:
            Translated from IDL: RJM - 2018/05/14
    """
    w_max     = numpy.max(w_km_vals[start:start + n_used])
    nw_fwd    = int(numpy.ceil(w_max / (2.0 * dx_km)))
    start_fwd = int(start + nw_fwd)
    n_fwd     = int(n_used - 2 * nw_fwd)

    # Compute product of wavenumber and RIP distance.
    kD_vals = 2.0 * numpy.pi * D_km_vals / lambda_km_vals

    # Compute Cosine of opening angle.
    cosb = numpy.cos(B_rad_vals)

    # Precompute variables for speed.
    cosbD = cosb/D_km_vals
    cosb2 = cosb*cosb

    # Compute cosine and sine of ring azimuth angle.
    cosphi0 = numpy.cos(phi_rad_vals)
    sinphi0 = numpy.sin(phi_rad_vals)

    # Precompute squares of variables for speed.
    dsq = D_km_vals*D_km_vals
    rsq = rho_km_vals*rho_km_vals

    # Define window function.
    fw = windows.FUNCTIONS[wtype]["func"]

    # Define normalization function.
    nrm = normalize_window.normalize_window

    # Set inverse function to FFT or Integration.
    if fft:
        finv = fresnel_inverse_fft
    else:
        finv = fresnel_inverse

    # Set psi approximation.
    if psitype == 'mtr4':
        psif = psi_quartic
    else:
        psif = psi_func

    # Create an empty array for reconstructed complex transmittance.
    T_hat_fwd_vals = 0.0 * T_vals

    # Set the first computed point to the 'start' variable.
    center = start_fwd

    # Compute first window width and window function.
    w_init = w_km_vals[center]
    w_func = fw(w_init, dx_km, error_check=False)

    # Compute number of points in window function
    nw = numpy.size(w_func)

    # Computed range about the first point
    crange = numpy.arange(int(center-(nw-1)/2), int(1+center+(nw-1)/2))-1

    # Compute first approximation for stationary phase.
    phi_s_rad = phi_rad_vals[center]

    # Compute Cosines and Sines of first approximation.
    cp = numpy.cos(phi_s_rad)
    sp = numpy.sin(phi_s_rad)

    # Factor used for first Newton-Raphson iteration
    dphi_fac = cosb2*cosphi0*sinphi0 / (1.0-cosb2*sinphi0*sinphi0)

    if psitype == 'fresnel':
        for i in numpy.arange(n_fwd):

            # Current point being computed.
            center = start_fwd+i

            # Ring radius of current point.
            r0 = rho_km_vals[center]

            # Window width for current point.
            w = w_km_vals[center]

            # Window function for current point.
            w_func = fw(w, dx_km, error_check=False)

            # Number of points in current window.
            nw = numpy.size(w_func)

            # Computed range of points.
            crange = numpy.arange(
                int(center-(nw-1)/2), int(1+center+(nw-1)/2)
            ) - 1

            # Computed ring radius range and Fresnel scale.
            r = rho_km_vals[crange]
            F = F_km_vals[center]

            # Compute psi for with stationary phase value
            psi_vals = (numpy.pi/2.0)*((r-r0)/F)*((r-r0)/F)

            # Compute kernel function for Fresnel inverse
            ker = w_func*numpy.exp(-1j*psi_vals)

            # Range of diffracted data that falls inside the window
            T = T_vals[crange]

            # Compute 'approximate' Fresnel Inversion for current point
            T_hat_fwd_vals[center] = finv(T_hat, ker, dx_km, F, error_check=False)

            # If normalization has been set, normalize the reconstruction
            if norm:
                T_hat_fwd_vals[center] *= nrm(dx_km, ker, F)

    else:
        for i in numpy.arange(n_fwd):

            # Current point being computed.
            center = start_fwd+i

            # Compute current radius and RIP distance.
            r0 = rho_km_vals[center]
            d2 = dsq[center]

            # Compute square of ring radius.
            r02 = rsq[center]
            phi0 = phi_rad_vals[center]

            # Precomputed variables for speed.
            cbd = cosbD[center]
            cp0 = cosphi0[center]
            sp0 = sinphi0[center]

            # Compute product of wavenumber and RIP Distance.
            kD = kD_vals[center]

            # Current window width and Fresnel scale.
            d = D_km_vals[center]
            b = B_rad_vals[center]
            w = w_km_vals[center]
            F = F_km_vals[center]

            # If the window width has changed, recompute variables.
            if numpy.abs(w_init - w) >= 2.0 * dx_km:

                # Reset w_init and recompute window function.
                w_init = w
                w_func = fw(w, dx_km, error_check=False)

                # Reset number of window points
                nw = numpy.size(w_func)

                # Computed range for current point
                crange = numpy.arange(
                    int(center-(nw-1)/2), int(1+center+(nw-1)/2)
                )

                # Ajdust ring radius by dx_km.
                r = rho_km_vals[crange]

                # Compute square of ring radius.
                r2 = rsq[crange]

                # First iteration of Newton-Raphson.
                dphi_s_rad = dphi_fac[center] * (r - r0) / r0

                # Compute new angle.
                phi_s_rad = phi_rad_vals[center] - dphi_s_rad

                # Compute Sines and Cosines of new angle.
                cp = numpy.cos(phi_s_rad)
                sp = numpy.sin(phi_s_rad)

            # If the window width has not changed, perform Newton-Raphson.
            else:

                # Adjust computed range by dx_km.
                crange += 1

                # Ajdust ring radius by dx_km.
                r = rho_km_vals[crange]

                # Compute square of ring radius.
                r2 = rsq[crange]

                # Compute Xi variable (MTR86 Equation 4b).
                xi = cbd * (r0 * cp0 - r * cp)

                # Compute Eta variable (MTR86 Equation 4c).
                eta = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

                # Compute intermediate variables for partial derivatives.
                v1 = r * cbd * sp
                v2 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v4 = numpy.sqrt(1.0 + 2.0*xi + eta)
                v5 = cbd * r * cp
                v6 = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2

                # Compute variables used for second partial derivative.
                dphia = (2.0*v5 + v6)/(2.0 * v4)
                dphib = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))

                # Compute First and Second Partial Derivatives of psi
                psi_d1 = (2.0*v1 + v2) / (2.0 * v4) - v1
                psi_d2 = dphia - dphib

                # Compute Newton-Raphson perturbation
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad

                # Compute Cosines and Sines new angle
                cp = numpy.cos(phi_s_rad)
                sp = numpy.sin(phi_s_rad)

            loop = 0

            # Perform Newton-Raphson on phi.
            while numpy.max(numpy.abs(dphi_s_rad)) > 1.E-10:

                # Compute Xi variable (MTR86 Equation 4b).
                xi = cbd * (r0 * cp0 - r * cp)

                # Compute Eta variable (MTR86 Equation 4c).
                eta = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2

                # Compute intermediate variables for partial derivatives.
                v1 = r * cbd * sp
                v2 = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / d2
                v4 = numpy.sqrt(1.0 + 2.0*xi + eta)
                v5 = cbd * r * cp
                v6 = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2

                # Compute variables used for second partial derivative.
                dphia = (2.0*v5 + v6)/(2.0 * v4)
                dphib = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))

                # Compute First and Second Partial Derivatives of psi
                psi_d1 = (2.0*v1 + v2) / (2.0 * v4) - v1
                psi_d2 = dphia - dphib

                # Compute Newton-Raphson perturbation
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad

                # Compute Cosines and Sines new angle
                cp = numpy.cos(phi_s_rad)
                sp = numpy.sin(phi_s_rad)

                # Add one to loop variable for each iteration
                loop += 1

                if loop > 5:
                    break

            # Compute psi for with stationary phase value
            psi_vals = kD * psif(r, r0, phi_s_rad, phi0, d, b, error_check=False)

            # Compute kernel function for Fresnel inverse
            ker = w_func*numpy.exp(1j*psi_vals)

            # Range of diffracted data that falls inside the window
            T = T_vals[crange]

            # Compute 'approximate' Fresnel Inversion for current point
            T_hat_fwd_vals[center] = finv(T, ker, dx_km, F, error_check=False)

            # If normalization has been set, normalize the reconstruction
            if norm:
                T_hat_fwd_vals[center] *= nrm(dx_km, ker, F)

            if verbose:
                print(f'Pt: {i} Tot: {n_fwd} Width: {nw} Psi Iters: {loop}')

    return T_hat_fwd_vals
