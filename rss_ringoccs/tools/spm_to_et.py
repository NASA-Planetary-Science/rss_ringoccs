#!/usr/bin/env python
"""

spm_to_et.py

Purpose: Calculate ephemeris time given a set of SPM values and appropriate
         kernels. Called inside of calc_f_sky_recon.py

Revisions:
      gjs_spm_to_et.py
   2018 Feb 22 - gsteranka - Original version
      spm_to_et.py
   2018 Mar 20 - gsteranka - Copy to official version and remove
                             debug steps
    2018 Jun 18 - jfong - remove 'from spiceypy' from import spicepy as spice
"""

import numpy as np
import spiceypy as spice

def spm_to_et(spm, doy, year, kernels=None):
    """
    Function to calculate ephemeris time from given SPM
    """

    # Leap seconds kernel
    if kernels is None:
        kernels = '../../../kernels/naif0012.tls'

    spice.kclear()
    spice.furnsh(kernels)

    n_pts = len(spm)

    et_sec_vals = np.zeros(n_pts)

    hours = (spm / 3600.0).astype(int)
    remainder_spm = (spm - hours*3600.0).astype(int)
    minutes = (remainder_spm / 60.0).astype(int)
    seconds = spm - hours*3600.0 - minutes*60.0

    # Test if it's a leap year
    if (year % 4) == 0:
        days_per_year = 366
    else:
        days_per_year = 365

    for i in range(n_pts):

        # Adjust if hour makes it go on to next day
        doy_adjusted = doy + hours[i]/24

        # Adjust if day goes on to next year
        this_year = (year + np.floor((doy_adjusted - 1)/days_per_year)).astype(int)
        this_doy = (1 + (doy_adjusted-1) % days_per_year).astype(int)

        utc_str = (str(this_year) + '-' + str(this_doy).zfill(3) + 'T'
            + str(hours[i] % 24).zfill(2) + ':'
            + str(minutes[i]).zfill(2) + ':'
            + '{:2.13f}'.format(seconds[i]))

        et_sec_vals[i] = spice.utc2et(utc_str)

    return et_sec_vals


if __name__ == '__main__':
    pass
