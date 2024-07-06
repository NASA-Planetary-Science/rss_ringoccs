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
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
# Pylint warns about numpy not having "reverse". This is false.
# pylint: disable = no-member

# Pylint also wants certain variables to be all caps. The original TC2017
# coding style didn't follow this. Ignore this warning.
# pylint: disable = invalid-name
import sys
import time
import numpy
from scipy.special import lambertw
from scipy.special import iv
import spiceypy as spice
import pandas as pd
import matplotlib.pyplot as plt

t1 = time.time()
data = 'data.txt'
df = pd.read_csv(data)

# Restore Parameters

# Requested Resolution (km)
res_km = 1.0

# Range to compute inversion
rrange = numpy.array([87410.0, 87610.0])

# Restore Variables

# Ring intercept point radius.
rho_km_vals = numpy.array(df.rho_km_vals)

# Uncorrected phase, radians.
phase_rad_vals = numpy.array(df.phase_rad_vals)

# Uncorrected power.
raw_power_vals = numpy.array(df.p_norm_vals)

# Ring-azimuth angle, radians.
phi_rad_vals = numpy.array(df.phi_rad_vals)

# Ring-opening angle, radians.
b_rad_vals = numpy.array(df.b_rad_vals)

# Sky frequency values, in Hz.
f_sky_hz_vals = numpy.array(df.f_sky_hz_vals)

# Spacecraft-rip distance, km.
d_km_vals = numpy.array(df.d_km_vals)

# drho/dt, in km/s
rho_dot_kms_vals = numpy.array(df.rho_dot_kms_vals)

# Distance between points, km.
dx_km = rho_km_vals[1]-rho_km_vals[0]
drho_dt = [numpy.min(rho_dot_kms_vals), numpy.max(rho_dot_kms_vals)]

# Number of points.
n = len(rho_km_vals)

# Normalized Equivalent Width.
neq = 1.651920

# USO Allen Deviation
sigma = 2.E-13

# Angular frequency (Hz)
omega = 2.0 * numpy.pi * f_sky_hz_vals

# Number of iterations for psi.
n_iter = 6

# Phase increment
dphi_rad = 0.00001
dphisq = dphi_rad * dphi_rad

# If dx is negative, then occultation is ingress and variables must reverse.
if dx_km < 0:
    rho_km_vals.reverse()
    phase_rad_vals.reverse()
    raw_power_vals.reverse()
    phi_rad_vals.reverse()
    b_rad_vals.reverse()
    f_sky_hz_vals.reverse()
    d_km_vals.reverse()

    # dr/dt must be positive
    rho_dot_kms_vals = -numpy.flipud(rho_dot_kms_vals)
    dx_km = -dx_km

# Compute other variables that are not contained in the original save file.

# Wavelength of signal.
lambda_sky = spice.clight() / f_sky_hz_vals

# Sine of ring opening angle.
mu = numpy.sin(numpy.abs(b_rad_vals))

# The Fresnel Scale
f_km_vals = numpy.sqrt(
    lambda_sky * d_km_vals * 0.5 * (
        1.0 - (numpy.square(numpy.cos(b_rad_vals) * numpy.sin(phi_rad_vals)))
    ) / numpy.square(numpy.sin(b_rad_vals))
)
T_hat_vals = numpy.sqrt(raw_power_vals) * numpy.exp(1j * phase_rad_vals)

# Create variable to determine if an error has occurred.
error_code = 0

# Make sure the requested resolution has double valued precision.
if not isinstance(res_km, float):
    print('Illegal res_km Type. Must be double valued.')
    error_code  +=  2

# Check to make sure there is no ambiguity about the type of occultation.
if drho_dt[0] < 0 < dx_km:
    print("rho_km_vals are increasing but rho_dot_kms_vals is negative")
    print("I can't determine if this is an ingress or an egress occultation")
    print("Correct so rho_dot_kms_vals and rho_km_vals[1]-rho_km_vals[0]")
    print("Have the same sign.")
    error_code += 4

# Check to make sure there is no ambiguity about the type of occultation.
if dx_km < 0 < drho_dt[1]:
    print("rho_km_vals are decreasing but rho_dot_kms_vals is positive")
    print("I can't determine if this is an ingress or an egress occultation")
    print("Correct so that rho_dot_kms_vals and rho_km_vals[1]-rho_km_vals[0]")
    print("Have the same sign.")
    error_code += 8

# Check that the number of phase points equals the number of radius points.
if len(phase_rad_vals) != n:
    print('Bad Input: len(phase_rad_vals) != len(rho_km_vals)')
    error_code +=  16

# Check that the number of power points equals the number of radius points.
if len(raw_power_vals) != n:
    print('Bad Input: len(raw_power_vals) != len(rho_km_vals)')
    error_code +=  32

# Check the the ring azimuth and ring radius have the same number of points.
if len(phi_rad_vals) != n:
    print('Bad Input: len(phi_rad_vals) != len(rho_km_vals)')
    error_code +=   64

# Check that the ring-opening and ring radius have the same number of points.
if len(b_rad_vals) != n:
    print('Bad Input: n_elements(B_rad_vals) ne n_elements(rho_km_vals)')
    error_code += 128

# Check that the sky frequency and ring radius have the same number of points.
if len(f_sky_hz_vals) != n:
    print('Bad Input: n_elements(f_sky_Hz_vals) ne n_elements(rho_km_vals)')
    error_code += 256

# Check the the D-values and the ring radius have the same number of points.
if len(d_km_vals) != n:
    print('Bad Input: n_elements(D_km_vals) ne n_elements(rho_km_vals)')
    error_code += 512

# Check rho_dot_km_vals and the ring radius have the same number of points.
if len(rho_dot_kms_vals) != n:
    print('Bad Input: n_elements(D_km_vals) ne n_elements(rho_km_vals)')
    error_code += 1024

# Make sure that ring radius is a double precision array.
if not isinstance(rho_km_vals, numpy.ndarray):
    print('Bad Input: rho_km_vals is not double')
    error_code += 2048

# Make sure that phase is a double precision array.
if not isinstance(phase_rad_vals, numpy.ndarray):
    print('Bad Input: phase_rad_vals is not double')
    error_code += 4096

# Make sure that phase is a double precision array.
if not isinstance(raw_power_vals, numpy.ndarray):
    print('Bad Input: raw_power_vals is not double')
    error_code += 8192

# Make sure the ring-azimuth is a double precision array
if not isinstance(phi_rad_vals, numpy.ndarray):
    print('Bad Input: phi_rad_vals is not double')
    error_code += 16384

# Make sure the B-values are a double precision array.
if not isinstance(b_rad_vals, numpy.ndarray):
    print('Bad Input: b_rad_vals is not double')
    error_code += 32678

# Make sure the sky-frequency is a double precision array.
if not isinstance(f_sky_hz_vals, numpy.ndarray):
    print('Bad Input: f_sky_hz_vals is not double')
    error_code += 65536

# Make sure the D-values are a double precision array.
if not isinstance(d_km_vals, numpy.ndarray):
    print('Bad Input: d_km_vals is not double')
    error_code += 131072

# Make sure the rho_dot_kms_vals are a double precision array.
if not isinstance(rho_dot_kms_vals, numpy.ndarray):
    print('Bad Input: rho_dot_kms_vals is not double')
    error_code += 262144

# If a problem did occur, stop the code and find the error.
if error_code != 0:
    sys.exit()

# Parameters for inverse of Resolution = f(Window)
alpha = 0.5 * numpy.square(omega * sigma) / rho_dot_kms_vals
P = res_km / (alpha * numpy.square(f_km_vals))

# The inverse exists only if P>1.
if min(P) < 1.0001:
    print('WARNING: Bad Points!')
    print(numpy.where(P<1.0001))
    print('Either rho_dot_kms_vals, F_km_vals, or res_km is to small.')
    sys.exit()

# Create window variable, window width (in km) for each point.
w_km_vals = numpy.abs(
    numpy.sqrt(2.0) * neq * (
        (P - 1.0)*lambertw(numpy.exp(P/(1-P))*P/(1-P))+P
    ) / (P - 1.0) / alpha
)

# Largest window used, km
w_max = numpy.abs(numpy.max(w_km_vals))

# If the window width is large, inversion will take a while.
if w_max > 1000:
    print('WARNING: Windows larger than 1000 km needed.')
    print('Max Window Width Used: ', w_max)

# Smallest possible starting point.
rho_min_lim = numpy.min(rho_km_vals)+numpy.ceil(w_max/2)

# Largest possible ending point.
rho_max_lim = numpy.max(rho_km_vals)-numpy.ceil(w_max/2)

# Beggining of requested range
rho_start = rho_km_vals[numpy.min(numpy.where(rho_km_vals > rrange[0]))]

# End of requested range
rho_end = rho_km_vals[numpy.max(numpy.where(rho_km_vals < rrange[1]))]

# Starting radius, taking constraints into consideration.
rho_min = numpy.max([rho_min_lim,rho_start])

# Ending radius, taking constraints into consideration.
rho_max = numpy.min([rho_max_lim,rho_end])

# Starting point.
start = numpy.min(numpy.where(rho_km_vals > rho_min))

# Ending point.
finish = numpy.max(numpy.where(rho_km_vals < rho_max))

# Actual range used.
range_km_vals = rho_km_vals[start:finish]

# Number of points used.
n_used = len(range_km_vals)

# Empty array for the corrected complex amplitude.
T_vals = numpy.zeros(n)*(1+1j)

# For matrix algebra
uni_m = numpy.array([[1],[1],[1],[1],[1]])

# Create variable for Lagrange Interpolation
dphi_val = numpy.array([[-2.0], [-1.0], [0.0], [1.0], [2.0]]) * dphi_rad

# The denominator of the window function.
window_factor = 1.0 / iv(0, 2.5 * numpy.pi)

# Calculate the corrected complex amplitude, point by point
for i in range(n_used):
    center = start + i
    n_pts_w = int(2.0 * numpy.ceil(w_km_vals[center] / (2.0*dx_km)) - 1.0)
    x = (numpy.array(range(n_pts_w)) - n_pts_w) * dx_km * 0.5
    window_arg = 2.5 * numpy.pi * numpy.sqrt(
        (1.0 - numpy.square(2.0 * x/w_km_vals[center]))
    )

    w_func = iv(0, window_arg) * window_factor
    nw = len(w_func)
    crange = numpy.array(range(int(center-(nw-1)/2),int(1+center+(nw-1)/2)))
    r = numpy.array([rho_km_vals[crange]])
    r0 = rho_km_vals[center]
    F_km_sq = f_km_vals[center] * f_km_vals[center]
    D = d_km_vals[center]
    B = b_rad_vals[center]
    f_sky_hz = f_sky_hz_vals[center]
    phi0 = phi_rad_vals[center]
    T_hat_Lu = numpy.array(T_hat_vals[crange])
    lambda_km = spice.clight()/f_sky_hz
    kD = (spice.twopi()/lambda_km)*D
    cos_phi0 = numpy.cos(phi0)
    sin_phi0 = numpy.sin(phi0)
    sin_phi0_sq = sin_phi0 * sin_phi0
    r0cos_phi0 = r0*cos_phi0
    Dsq = D*D
    r0_sq = r0*r0
    cosB  = numpy.cos(B)
    cosB_sq = cosB*cosB
    factor = (((cosB_sq)*cos_phi0*sin_phi0) / (1.0-(cosB_sq)*sin_phi0_sq)) / r0
    idd  = numpy.array([numpy.zeros(nw) + 1])
    r_id = numpy.dot(uni_m, r)
    r_sq = r_id * r_id
    phi_s_rad = numpy.array(phi0-factor*(r-r0))
    dphi = numpy.dot(dphi_val, idd)

    for kk in range(n_iter):
        phi_s_rad_vals = numpy.dot(uni_m,phi_s_rad) + dphi
        xi     = cosB*(r0cos_phi0-(r_id*numpy.cos(phi_s_rad_vals)))/D
        eta    = (r0_sq+r_sq-2.*r0*(r_id*numpy.cos(phi_s_rad_vals-phi0)))/Dsq
        FF     = kD*(numpy.sqrt(1.+2.*xi+eta)-(1.+xi))
        Psi_x  = (FF[0,::]-8.*FF[1,::]+8.*FF[3,::]-FF[4,::])/(12.*dphi_rad)
        Psi_xx = (-FF[0,::]+16.*FF[1,::]-30.*FF[2,::]+16.*FF[3,::]-FF[4,::])/\
        (12.*dphisq)
        dphi_s_rad = -Psi_x/Psi_xx
        phi_s_rad += dphi_s_rad

    xi = (cosB/D) * (r0*cos_phi0-r*numpy.cos(phi_s_rad))
    eta = (r0 * r0+ r*r - 2.0*r*r0*numpy.cos(phi_s_rad-phi0)) / Dsq
    psi_vals = numpy.reshape(kD*(numpy.sqrt(1.+2.*xi+eta)-(1.+xi)), nw)
    arg = numpy.exp(1j * psi_vals)
    F_km = f_km_vals[center]
    conv_kernel = arg*w_func
    fft_t_hat = numpy.fft.fftshift(numpy.fft.fft(T_hat_Lu))
    fft_conv = numpy.fft.fftshift(numpy.fft.fft(conv_kernel))
    inv_t_hat = numpy.fft.ifftshift(numpy.fft.ifft(fft_t_hat*fft_conv))
    inv_t_hat *= dx_km*(numpy.complex(1.,1.))/(2.*F_km)
    T_vals[center] = inv_t_hat[int((nw-1) / 2)]
    print(i, n_used - 1)

r = rho_km_vals[start:start + n_used]
power = numpy.square(numpy.abs(T_vals[start:start + n_used]))
tau = -mu[start:start + n_used] * numpy.log(power)
plt.plot(r, power)
t2 = time.time()

print('Hours:', int(numpy.floor((t2-t1)/3600)))
print('Minutes:',int(numpy.floor(t2-t1)/60))
print('Seconds:', int(numpy.mod(t2-t1,60)))
