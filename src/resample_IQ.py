#!/usr/bin/env python
"""

resample_IQ.py

Purpose: Resample I and Q from uniformly spaced time to uniformly spaced
         radius. Normally, you resample from raw resolution data.

Revisions:
      gjs_resample_IQ.py
   Feb 27 2018 - gsteranka - Original version
   Mar 01 2018 - gsteranka - Took pre_resample function outside of resample_IQ
                             function, since I found out it doesn't need to be
                             inside to be recognized in a call to resample_IQ
                             from another file
   Mar 07 2018 - gsteranka - Input is a complex number IQ_c, instead of
                             separate real and imaginary parts
   Mar 08 2018 - gsteranka - Change input arrays so rho starts within dr_km_tol
                             of an integer number of dr_km. Makes it so after
                             resampling, everything is at a more even number
      resample_IQ.py
   Mar 20 2018 - gsteranka - Copy to official version and remove debug steps
   Apr 09 2018 - gsteranka - Edited pre_resample function to handle ingerss
                             files as well as egress
"""

import numpy as np
from scipy import signal
from scipy.interpolate import interp1d

def pre_resample(rho_km, vec, freq):
    """Sub-function inside of resample_IQ to put vectors to uniform spaced
    radius at a similar spacing as raw resolution.  For ingress occultations,
    this step implicitly reverses the radius scale when interpolating"""

    # Average radius spacing over region
    ts_avg = abs(rho_km[-1]  - rho_km[0])/float(len(rho_km) - 1)
    p = 1
    q = int(round(1.0/(ts_avg*freq)))
    dr_grid = float(p)/(q*freq)

    # Uniform radius grid at near-raw resolution to interpolate to
    # For ingress, this implicitly reverses radius scale!
    n_pts = round(abs(rho_km[-1] - rho_km[0])/dr_grid)
    rho_grid = min(rho_km) + dr_grid*np.arange(n_pts)

    # Interpolate to near-raw resolution
    # For ingress, this implicitly reverses radius scale!
    vec_grid_interp = interp1d(rho_km, vec, kind='linear',
                               fill_value='extrapolate')
    vec_grid = vec_grid_interp(rho_grid)

    return rho_grid, vec_grid, p, q


def resample_IQ(rho_km, IQ_c, dr_desired, dr_km_tol=0.01, TEST=False):
    """Function to resample I and Q to uniformly spaced radius. Based off of
    Matlab's resample function

    Example:
        >>> fit_inst = FreqOffsetFit(RSR_inst, geo_inst, f_spm, f_offset, \
                f_USO, kernels, k=k, rho_exclude=rho_exclude)
        >>> (spm_raw, IQ_c_raw) = fit_inst.get_IQ_c()
        >>> (rho_km_desired, IQ_c_resampled) = resample_IQ(rho_km_raw, \
                IQ_c_raw, dr_desired, dr_km_tol=dr_km_tol, TEST=TEST)

    Args:
        rho_km (np.ndarray): Set of ring intercept point values at initial
            resolution before resampling
        IQ_c (np.ndarray): Frequency corrected complex signal at initial
            resolution before resampling
        dr_desired (float): Desired final radial spacing to resample to
        dr_km_tol (float): Maximum distance from an integer number of
            dr_desired that the first rho value will be at. Makes the final
            set of rho values more regular-looking. For example, if you say
            dr_km_tol=0.01 with a dr_desired of 0.25, your final set of rho
            values might look something like [70000.26, 70000.51, ...]
        TEST (bool): Testing variable to print out the first few resampled
            results
    """

    # Begin raw resolution rho within dr_km_tol of integer number of dr_desired
    rho_km_remainder = rho_km % dr_desired
    rho_km_start_ind = int((np.argwhere(rho_km_remainder < dr_km_tol))[0])
    rho_km = rho_km[rho_km_start_ind:-1]
    IQ_c = IQ_c[rho_km_start_ind:-1]

    I_c = np.real(IQ_c)
    Q_c = np.imag(IQ_c)

    # Pre-resampling steps. Interpolates to uniform radius at near-raw spacing
    (rho_km_uniform, I_c_uniform, p, q) = pre_resample(rho_km, I_c,
                                                       1.0/dr_desired)
    (rho_km_uniform, Q_c_uniform, p, q) = pre_resample(rho_km, Q_c,
                                                       1.0/dr_desired)

    # Downsample by factor q to desired final spacing
    I_c_desired = signal.resample_poly(I_c_uniform, p, q)
    Q_c_desired = signal.resample_poly(Q_c_uniform, p, q)

    rho_km_desired = rho_km_uniform[0] + dr_desired*np.arange(len(I_c_desired))

    if TEST:
        print('First 10 rho, I_c, Q_c:')
        for i in range(10):
            print('%24.16f %24.16f %24.16f' %
                  (rho_km_desired[i], I_c_desired[i], Q_c_desired[i]))

    return rho_km_desired, I_c_desired + 1j*Q_c_desired


if __name__ == '__main__':
    pass
