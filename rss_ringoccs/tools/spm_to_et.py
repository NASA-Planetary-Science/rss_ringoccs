"""

:Purpose:
    Calculate ephemeris time given a set of SPM values and appropriate kernels.

:Dependencies:
    #. numpy
    #. spiceypy
    #. sys
"""
import numpy as np
import spiceypy as spice
import sys

def spm_to_et(spm, doy, year, kernels=None):
    """
    Convert seconds past midnight to ephemeris seconds past J2000.

    Arguments
        :spm (*np.ndarray*): SPM values
        :doy (*int*): Day of year of observation
        :year (*int*): Year of observation

    Keyword Arguments
        :kernels (*str*): String specifying the appropriate ephemeris
                kernel file. If ``None``, sets the kernel file
                to
                ../../kernels/naif/CASSINI/kernels/lsk/naif0012.tls
                Default is ``None``.
    """

    if doy == 0:
        sys.exit('WARNING (spm_to_et): input doy is 0!')

    if year == 0:
        sys.exit('WARNING (spm_to_et): input year is 0!')

    year = int(year)
    doy = int(doy)


    # Leap seconds kernel
    if kernels is None:
        kernels = '../../kernels/naif/CASSINI/kernels/lsk/naif0012.tls'

    spice.kclear()
    spice.furnsh(kernels)

    if isinstance(spm, float):
        hours = int(spm/3600.)
        remainder_spm = int(spm - hours*3600.0)
        minutes = int(remainder_spm/60.)
        seconds = spm - hours*3600.0 - minutes*60.0
        n_pts = 1
        hours = [hours]
        minutes = [minutes]
        seconds = [seconds]
    else:
        hours = (spm / 3600.0).astype(int)
        remainder_spm = (spm - hours*3600.0).astype(int)
        minutes = (remainder_spm / 60.0).astype(int)
        seconds = spm - hours*3600.0 - minutes*60.0
        n_pts = len(spm)


    et_sec_vals = np.zeros(n_pts)

    # Check for leap year
    if (year % 4) == 0:
        days_per_year = 366
    else:
        days_per_year = 365

    for i in range(n_pts):

        # Adjust if hour makes it go on to next day
        doy_adjusted = doy + hours[i]/24

        # Adjust if day goes on to next year
        if n_pts == 1:
            this_year = int(year +
                np.floor((doy_adjusted - 1)/days_per_year))
            this_doy = int(1 + (doy_adjusted-1) % days_per_year)
        else:
            this_year = (year +
                np.floor((doy_adjusted - 1)/days_per_year)).astype(int)
            this_doy = (1 + (doy_adjusted-1) % days_per_year).astype(int)

        utc_str = (str(this_year) + '-' + str(this_doy).zfill(3) + 'T'
            + str(hours[i] % 24).zfill(2) + ':'
            + str(minutes[i]).zfill(2) + ':'
            + '{:2.13f}'.format(seconds[i]))

        et_sec_vals[i] = spice.utc2et(utc_str)

    return et_sec_vals

"""
Revisions:
"""
