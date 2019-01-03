"""
Purpose:
    Calculate ephemeris time given a set of SPM values and appropriate
    kernels. Called by ``calc_f_sky_recon.py``.
"""

import numpy as np
import spiceypy as spice
import sys
import pdb

def spm_to_et(spm, doy, year, kernels=None):
    """
    Arguments:
        :spm (*np.ndarray*): SPM values
        :doy (*int*): Day of year of observation
        :year (*int*): Year of observation

    Keyword Arguments:
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




    # Leap seconds kernel
    if kernels is None:
        kernels = '../../kernels/naif/CASSINI/kernels/lsk/naif0012.tls'

    spice.kclear()
    spice.furnsh(kernels)

    hours = (spm / 3600.0).astype(int)
    remainder_spm = (spm - hours*3600.0).astype(int)
    minutes = (remainder_spm / 60.0).astype(int)
    seconds = spm - hours*3600.0 - minutes*60.0

    try:
        n_pts = len(spm)
    except TypeError:
        n_pts = len([spm])
        hours = [hours]
        minutes = [minutes]
        seconds = [seconds]

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
      gjs_spm_to_et.py
   2018 Feb 22 - gsteranka - Original version
      spm_to_et.py
   2018 Mar 20 - gsteranka - Copy to official version and remove
                             debug steps
    2018 Jul 23 - jfong - add input error checks (doy and year can't be 0)
    2018 Jul 30 - jfong - allow float inputs for spm (in addition to array)
"""
