#!/usr/bin/env python
"""

:Purpose:
    Resample I and Q from uniformly spaced time to uniformly spaced
    radius. This is set up to downsample from the raw resolution data.
"""

import numpy as np
from scipy import signal
from scipy.interpolate import interp1d


def pre_resample(rho_km, vec, freq):
    """
    Set vector sampling to be uniform with respect to radius at
    a spacing comparable to that of raw resolution.  For ingress
    occultations, this step implicitly reverses the radius scale
    when interpolating.

    Arguments
        :rho_km (*np.ndarray*): radius in kilometers
        :vec (*np.ndarray*): a single vector component I or Q of the
                complex signal
        :freq (*float*): radial sampling frequency

    Returns
        :rho_grid (*np.ndarray*): Radii at uniform spacing at which the
                 signal component is resampled
        :vec_grid (*np.ndarray*): Signal resampled with respect to radius
                 with a uniform spacing
        :p (*float*): upsampling rate to be used by scipy.signal.resample_poly.
                 This will always be unity because no upsampling is done.
        :q (*float*): downsampling rate to be used by
                  scipy.signal.resample_poly. This depends on the uniform
                  radial sampling rate at which `rho_grid` and `vec_grid` are
                  sampled.

    """
    # initial radius and radial range
    r0 = round(rho_km[0])
    dr = abs(rho_km[-1] - r0)

    # Average radius spacing over region
    ts_avg = dr / float(len(rho_km) - 1)
    # check to make sure desired resolution is allowed
    # if average radial sampling is larger than desired
    # radial resolution, then the desired resolution
    # is not achievable with this data set
    if ts_avg > 1./freq :
        # round up to nearest 0.025
        res = round(40*ts_avg)/40
        if res < ts_avg:
            res += 0.025
        # print warning + advisory for 16 kHz files
        print('\tWARNING:')
        warn_str = ('Desired DLP resolution not achievable.\n\t\t'
              +'Average raw radial sampling of '+str(round(ts_avg,3))
              +' km\n\t\tis larger than the '+str(1./freq)+' km value '
              +'assigned \n\t\tto the variable dr_km. Setting DLP '
              +'\n\t\tresolution to lowest achievable raw\n\t\tsampling of '
              +str(res)+' km.\n\t\t')
        advc_str = ('To achieve the desired '+str(1./freq)+' km resolution,\n\t\t'
              +'consider using the 16 khz data file for this\n\t\t'
              +'occultation.')
        print('\t\t'+warn_str+advc_str)
        # update spatial frequency
        freq = 1./res
    else:
        res = round(1./freq,3)
    # upsampling factor -- default is 1
    p = 1
    # downsampling factor
    q = int(round(1.0 / (ts_avg * freq)))
    #print(freq,1/freq,dr,ts_avg,p,q)
    dr_grid = float(p) / (float(q) * freq)

    # Uniform radius grid at near-raw resolution to which to interpolate.
    #     For ingress, this implicitly reverses radius scale!
    n_pts = round(dr / dr_grid)
    rho_grid = r0 + dr_grid * np.arange(n_pts)

    # Interpolate to near-raw resolution. For ingress, this implicitly
    #     reverses radius scale!
    vec_grid_interp = interp1d(rho_km, vec, kind='linear',
        fill_value='extrapolate')
    vec_grid = vec_grid_interp(rho_grid)

    return rho_grid, vec_grid, p, q, res


def resample_IQ(rho_km, IQ_c, dr_desired, verbose=False):
    """
    Resample I and Q to uniformly spaced radius. Based off of
    Matlab's ``resample`` function

    Arguments
        :rho_km (*np.ndarray*):
            Set of ring intercept point values at initial
            resolution before resampling
        :IQ_c (*np.ndarray*):
            Frequency-corrected complex signal at initial
            resolution before resampling
        :dr_desired (*float*):
            Desired final radial sample spacing
        :verbose (*bool*):
            Testing variable to print out the first few resampled
            results

    Returns
        :rho_km_desired (*np.ndarray*): array of ring radius at final desired
                spacing
        :IQ_c_desired (*np.ndarray*):
    """

    rho_km_diff = np.diff(rho_km)

    # If you didn't select indices right in a chord occultation, then
    #     rho_km might not be monotonically increasing or decreasing
    if np.any(rho_km_diff > 0) & np.any(rho_km_diff < 0):
        print('WARNING (resample_IQ.py): This routine assumes the rho_km')
        print('    input to be either monotonically increasing or')
        print('    monotonically decreasing. Current input has both')

    # Reverse ingress to be increasing radius. This lets the first radius
    #     be selected so that the final radii are at integer numbers of
    #     requested spacing
    if rho_km_diff[0] < 0:
        if verbose:
            print('DETECTED INGRESS (resample_IQ.py): reversing arrays')
        rho_km = rho_km[::-1]
        IQ_c = IQ_c[::-1]

    I_c = np.real(IQ_c)
    Q_c = np.imag(IQ_c)

    # Pre-resampling steps. Interpolates to uniform radius at near-raw spacing
    rho_km_uniform, I_c_uniform, p, q, dr = pre_resample(rho_km, I_c,
        1.0 / dr_desired)
    rho_km_uniform, Q_c_uniform, p, q, dr = pre_resample(rho_km, Q_c,
        1.0 / dr_desired)

    # Downsample by factor q to desired final spacing
    I_c_desired = signal.resample_poly(I_c_uniform, p, q)
    Q_c_desired = signal.resample_poly(Q_c_uniform, p, q)

    rho_km_desired = (rho_km_uniform[0]
        + dr * np.arange(len(I_c_desired)))

    IQ_c_desired = I_c_desired + 1j * Q_c_desired
    return rho_km_desired, IQ_c_desired

"""
Revisions:
"""
