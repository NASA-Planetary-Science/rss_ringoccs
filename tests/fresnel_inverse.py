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
#       Small experiment with the Fresnel inverse for a square well.           #
################################################################################
#   Author:     Ryan Maguire                                                   #
#   Date:       2025/04/24                                                     #
################################################################################
"""
import numpy

# Let's keep this constant for all of the modeling. 15 degrees, halfway
# through a Saturn year.
OPENING_ANGLE = 2.617993877991494E-01
COS_B = numpy.cos(OPENING_ANGLE)
SIN_B = numpy.sin(OPENING_ANGLE)
COS_B_SQUARED = COS_B * COS_B
SIN_B_SQUARED = SIN_B * SIN_B
RCPR_SIN_B_SQUARED = 1.0 / SIN_B_SQUARED

# For now, let's keep this constant too. In future scripts we can let it vary.
DISTANCE = 2.7E+05
RCPR_DISTANCE = 1.0 / DISTANCE
RCPR_DISTANCE_SQUARED = RCPR_DISTANCE * RCPR_DISTANCE

# Common factor used, ratio of cos(B) and the distance D.
ALPHA = COS_B * RCPR_DISTANCE

# Wavenumber for a radio wave with a 3.6 cm wavelength.
WAVENUMBER = 1.74532925199432957692E+05

# The scale factor for the Fresnel kernel is kD.
SCALE_FACTOR = WAVENUMBER * DISTANCE

# Maximum number of iterations allowed in Newton's method.
MAX_ITERS = 8

# The tolerance, or epsilon value, for Newton's method.
EPSILON = 1.0E-8

# Scale factor that appears in the Fresnel scale.
PI = 3.141592653589793
FRESNEL_FACTOR = PI * DISTANCE / (WAVENUMBER * RCPR_SIN_B_SQUARED)

# Requested resolution, in kilometers, for the forward models.
RESOLUTION = 1.0

# Computes the (negative of the) xi factor from MTR86.
def xi_factor(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            xi_factor
        Purpose:
            Computes the negative of the xi factor that appears in MTR86.
            We adopt the convention of negating xi since this makes the
            generating function for the Legendre polynomials appear more
            naturally when deriving the stationary azimuth approximation.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    cos_phi = numpy.cos(phi_vals)
    cos_phi0 = numpy.cos(phi0_vals)
    return ALPHA * (r_vals * cos_phi - r0_vals * cos_phi0)

# Computes the eta factor from MTR86.
def eta_factor(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            eta_factor
        Purpose:
            Computes the eta factor that appears in MTR86. Unlike the xi
            function above, the eta factor here agrees with the MTR86 paper.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    cos_diff = numpy.cos(phi_vals - phi0_vals)
    numerator = r_vals*r_vals + r0_vals*r0_vals - 2.0*r_vals*r0_vals*cos_diff
    return numerator * RCPR_DISTANCE_SQUARED

# Computes the partial derivative of xi with respect to phi.
def dxi_dphi(r_vals, phi_vals):
    """
        Function:
            dxi_dphi
        Purpose:
            Computes d xi / d phi.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
    """
    sin_phi = numpy.sin(phi_vals)
    return -ALPHA * r_vals * sin_phi

# Computes the partial derivative of eta with respect to phi.
def deta_dphi(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            deta_dphi
        Purpose:
            Computes d eta / d phi.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    sin_diff = numpy.sin(phi_vals - phi0_vals)
    return 2.0 * r0_vals * r_vals * sin_diff * RCPR_DISTANCE_SQUARED

# Computes the second partial derivative of xi with respect to phi.
def d2xi_dphi2(r_vals, phi_vals):
    """
        Function:
            d2xi_dphi2
        Purpose:
            Computes d^2 xi / d phi^2.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
    """
    cos_phi = numpy.cos(phi_vals)
    return -ALPHA * r_vals * cos_phi

# Computes the second partial derivative of eta with respect to phi.
def d2eta_dphi2(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            deta_dphi
        Purpose:
            Computes d^2 eta / d phi^2.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    cos_diff = numpy.cos(phi_vals - phi0_vals)
    return 2.0 * r0_vals * r_vals * cos_diff * RCPR_DISTANCE_SQUARED

# Computes psi, the Fresnel kernel.
def psi(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            psi
        Purpose:
            Computes the Fresnel kernel.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    xi_vals = xi_factor(r_vals, r0_vals, phi_vals, phi0_vals)
    eta_vals = eta_factor(r_vals, r0_vals, phi_vals, phi0_vals)
    main_term = numpy.sqrt(1.0 - 2.0 * xi_vals + eta_vals)
    return SCALE_FACTOR * (main_term - 1.0 + xi_vals)

# Computes the partial derivative of psi with respect to phi.
def dpsi_dphi(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            dpsi_dphi
        Purpose:
            Computes d psi / d phi
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    xi_vals = xi_factor(r_vals, r0_vals, phi_vals, phi0_vals)
    dxi_vals = dxi_dphi(r_vals, phi_vals)

    eta_vals = eta_factor(r_vals, r0_vals, phi_vals, phi0_vals)
    deta_vals = deta_dphi(r_vals, r0_vals, phi_vals, phi0_vals)

    numer = -2.0*dxi_vals + deta_vals
    denom = 2.0 * numpy.sqrt(1.0 - 2.0 * xi_vals + eta_vals)

    return SCALE_FACTOR * (numer / denom + dxi_vals)

# Computes the second partial derivative of psi with respect to phi.
def d2psi_dphi2(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            d2psi_dphi2
        Purpose:
            Computes d^2psi / d phi^2
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    xi_vals = xi_factor(r_vals, r0_vals, phi_vals, phi0_vals)
    dxi_vals = dxi_dphi(r_vals, phi_vals)
    d2xi_vals = d2xi_dphi2(r_vals, phi_vals)

    eta_vals = eta_factor(r_vals, r0_vals, phi_vals, phi0_vals)
    deta_vals = deta_dphi(r_vals, r0_vals, phi_vals, phi0_vals)
    d2eta_vals = d2eta_dphi2(r_vals, r0_vals, phi_vals, phi0_vals)

    main_term = numpy.sqrt(1.0 - 2.0 * xi_vals + eta_vals)
    main_term_cubed = main_term * main_term * main_term

    left = (-2.0 * d2xi_vals + d2eta_vals) / (2.0 * main_term)
    right = numpy.square(-2.0 * dxi_vals + deta_vals) / (4.0 * main_term_cubed)

    return SCALE_FACTOR * (left - right + d2xi_vals)

# Computes the stationary azimuth angle phi_s.
def stationary_angle(r_vals, r0_vals, phi_vals, phi0_vals):
    """
        Function:
            d2psi_dphi2
        Purpose:
            Computes the angle phi_s such that d psi / d phi = 0.
        Argument:
            r_vals (numpy.ndarray):
                The dummy variable, often integrated over.
            r0_vals (numpy.ndarray):
                The radius of the point of interest.
            phi_vals (numpy.ndarray):
                The azimuth angles corresponding to the r_vals array.
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """
    phi_s = phi_vals

    for _ in range(MAX_ITERS):
        dpsi_vals = dpsi_dphi(r_vals, r0_vals, phi_s, phi0_vals)
        d2psi_vals = d2psi_dphi2(r_vals, r0_vals, phi_s, phi0_vals)

        phi_s = phi_s - dpsi_vals / d2psi_vals

        if numpy.max(numpy.abs(dpsi_vals)) < EPSILON:
            break

    return phi_s

# Computes the Fresnel scale from the azimuth angle.
def fresnel_scale(phi0_vals):
    """
        Function:
            fresnel_scale
        Purpose:
            Computes F, the Fresnel scale.
        Argument:
            phi0_vals (numpy.ndarray):
                The angle for r0_vals.
    """

    sin_phi0 = numpy.sin(phi0_vals)
    sin_phi0_squared = sin_phi0 * sin_phi0
    f_squared = FRESNEL_FACTOR * (1.0 - COS_B_SQUARED * sin_phi0_squared)
    return numpy.sqrt(f_squared)

def window_width(resolution, phi0_vals):
    """
        Computes the window width needed for a given resolution.
        This takes the Fresnel scale into account.
    """
    f_scale = fresnel_scale(phi0_vals)
    return 2.0 * f_scale * f_scale / resolution

# Square well for diffraction modeling.
def square_well(left, right, x_vals):
    """
        Computes a simple square well.
    """
    left_edge = 1.0 - numpy.heaviside(x_vals - left, 0.5)
    right_edge = numpy.heaviside(x_vals - right, 0.5)
    return left_edge + right_edge
