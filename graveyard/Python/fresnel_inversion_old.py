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
#       Computes the Fresnel inversion of a diffracted data set.               #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""
# PDS names are all caps. Pylint doesn't like this.
# pylint: disable = invalid-name
import numpy
from . import normalize_window
from . import windows

def fresnel_inversion(rho_vals, F_vals, phi_rad_vals, B_vals,
                      D_vals, T_hat_vals, lambda_vals, w_vals,
                      dx, wtype, start, n_used, norm = True,
                      fft = False, verbose = True, psitype = 'full'):
    """
        Function:
            fresnel_inversion
        Purpose:
            Computes the fresnel inversion from a set of diffracted data
            using a 'fast' method to speed up computation time. This is
            achieved by computing cosine and sine function in the outer for
            loop, and then passing these computed values into the functions
            that need them. The normal version passes the arguments to the
            functions, and then cosines and sines are computed within the
            function. For small data sets or coarser resolutions, the normal
            version is faster. Both the normal and fast versions output
            completely identical results.
        Variables:
            rho_vals:
                Ring radius, in kilometers.
            F_vals:
                Fresnel scale, in kilometers.
            phi_rad_vals:
                Ring azimuth angle, in radians.
            B_rad_vals:
                Ring opening angle, in radians.
            lambda_sky_vals:
                Wavelength of recieved signal, in kilometers.
            D_vals:
                Spacecraft-RIP distance, in kilometers.
            dx:
                Sampling spacing, in kilometers.
            T_vals:
                Reconstructed complex transmittance.
            w_vals:
                Window width, in kilometers.
            wtype:
                Window used in reconstruction, string.
            start:
                Starting point of reconstructed data.
            n_used:
                Number of reconstructed points.
        Keywords:
            Normalize:
                Parameter for normalizing the complex transmittance by
                the window function that is used. Default is True. Set to False
                to skip this feature.
        Output:
            T_vals:
                Reconstructed Complex transmittance.
    """
    # Compute necessary variables.
    kD_vals = 2.0*numpy.pi*D_vals/lambda_vals

    # Define functions
    fw = windows.FUNCTIONS[wtype]["func"]
    nrm = normalize_window.normalize_window

    if fft:
        finv = fresnel_inverse_fft
    else:
        finv = fresnel_inverse

    # Calculate the corrected complex amplitude, point by point
    T_vals = T_hat_vals * 0.0
    w_init = w_vals[start]
    w_func = fw(w_init,dx)
    nw = numpy.size(w_func)
    crange = numpy.arange(int(start-(nw-1)/2),int(1+start+(nw-1)/2))
    phi_s_rad1 = phi_rad_vals[start]

    if psitype == 'taylor2':
        for i in numpy.arange(n_used):
            center      = start+i
            r0          = rho_vals[center]
            w           = w_vals[center]
            w_func      = fw(w,dx,error_check=False)
            nw          = numpy.size(w_func)
            crange      = numpy.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
            r           = rho_vals[crange]
            F           = F_vals[center]
            x           = (r-r0)/F
            psi_vals    = (numpy.pi/2.0)*x*x
            ker         = w_func*numpy.exp(-1j*psi_vals)
            T_hat       = T_hat_vals[crange]

            T_vals[center] = finv(T_hat,ker,dx,F)

            if norm:
                T_vals[center] *= nrm(r, w_func, F)

    elif psitype in ["mtr2", "mtr3", "mtr4"]:
        for i in numpy.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            d      = D_vals[center]
            b      = B_vals[center]
            phi0   = phi_rad_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            F      = F_vals[center]
            if numpy.abs(w_init - w)>= 2.0*dx:
                w_init     = w
                w_func     = fw(w,dx,error_check=False)
                nw         = numpy.size(w_func)
                crange     = numpy.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                dphi_s_rad = psi_factor(r,r0,b,phi0,error_check=False)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad

            else:
                crange    += 1
                r          = rho_vals[crange]
                phi_s_rad  = phi_s_rad1
                psi_d1     = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2     = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad

            loop = 0

            # Perform Newton-Raphson on phi.
            while numpy.max(numpy.abs(dphi_s_rad)) > 1.E-8:
                psi_d1     = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2     = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                loop      += 1
                if loop > 5:
                    break
            phi_s_rad1 = phi_s_rad

            # Compute psi and then compute averages and differences across psi.
            psi_full   = kD * psi(r,r0,d,b,phi_s_rad,phi0,error_check=True)
            n1 = 0                                      # Left Endpoint.
            n2 = numpy.min((r0-w/4<=r).nonzero())       # Left midpoint.
            n3 = numpy.max((r0+w/4>=r).nonzero())       # Right midpoint.
            n4 = nw-1                                   # Right endpoint.
            d_psi_half = psi_full[n3]-psi_full[n2]      # Midpoint difference.
            d_psi_full = psi_full[n4] - psi_full[n1]    # Endpoint difference.
            a_psi_half = (psi_full[n3]+psi_full[n2])/2  # Midpoint average.
            a_psi_full = (psi_full[n1]+psi_full[n4])/2  # Endpoint average.
            x = r - r0
            w = numpy.max(x) - numpy.min(x)

            #Compute coefficients for the polynomial expansion
            c1 = (8.0*d_psi_half-d_psi_full)/(3.0*w)            # Linear term
            c2 = 4.0*(16.0*a_psi_half-a_psi_full)/(3.0*w*w)     # Quadratic term
            c3 = 16.0*(d_psi_full-2.0*d_psi_half)/(3.0*w*w*w)   # Cubic term
            c4 = 64.0*(a_psi_full-4.0*a_psi_half)/(3.0*w*w*w*w) # Quartic term

            psi_vals = c1*x + c2*x*x        #Second order appoximation
            if psitype == 'mtr3':
                psi_vals += c3*x*x*x        #Third order approximation

            elif psitype == 'mtr4':
                psi_vals += (c3+c4*x)*x*x*x #Fourth order approximation

            ker = w_func*numpy.exp(-1j*psi_vals)
            T_hat = T_hat_vals[crange]

            T_vals[center] = finv(T_hat, ker, dx, F)

            if norm:
                T_vals[center] *= nrm(r, w_func, F)

            if verbose:
                print(f'Pt: {i} Tot: {n_fwd} Width: {nw} Psi Iters: {loop}')

    elif psitype == 'full':
        for i in numpy.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            d      = D_vals[center]
            b      = B_vals[center]
            phi0   = phi_rad_vals[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            F      = F_vals[center]
            if numpy.abs(w_init - w) >= 2.0*dx:
                w_init     = w
                w_func     = fw(w,dx,error_check=False)
                nw         = numpy.size(w_func)
                crange     = numpy.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                dphi_s_rad = psi_factor(r,r0,b,phi0,error_check=False)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
            else:
                crange    += 1
                r          = rho_vals[crange]
                phi_s_rad  = phi_s_rad1
                psi_d1     = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2     = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad

            loop = 0

            # Perform Newton-Raphson on phi.
            while numpy.max(numpy.abs(dphi_s_rad)) > 1.E-10:
                psi_d1 = psi_d1_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                psi_d2 = psi_d2_phi(r,r0,d,b,phi_s_rad,phi0,error_check=False)
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad

                loop += 1

                if loop > 5:
                    break

            phi_s_rad1 = phi_s_rad

            # Compute psi and then compute averages and differences across psi.
            psi_vals = kD*psi(r, r0, d, b, phi_s_rad, phi0)
            ker = w_func*numpy.exp(-1j*psi_vals)
            T_hat = T_hat_vals[crange]

            T_vals[center] = finv(T_hat,ker,dx,F)
            if norm:
                T_vals[center] *= nrm(r, psi_vals, w_func, F)

            if verbose:
                print(f'Pt: {i} Tot: {n_fwd} Width: {nw} Psi Iters: {loop}')

    else:
        raise TypeError(f"Illegal psitype: {psitype}")

    return T_vals
