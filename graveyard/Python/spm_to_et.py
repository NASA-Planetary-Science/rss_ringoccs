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
#       Calculate ephemeris time from SPM values and appropriate kernels.      #
################################################################################
#   Author: Glenn Steranka                                                     #
#   Date:   2018/02/22                                                         #
################################################################################
#                               Revision History                               #
################################################################################
#   2018/07/23: Jolene Fong                                                    #
#       Added input error checks (doy and year can't be 0).                    #
#   2018/07/30: Jolene Fong                                                    #
#       Allow float inputs for SPM.                                            #
################################################################################
"""
# Disabling a few pylint warnings.
import numpy
import spiceypy

def spm_to_et(spm, doy, year, kernels = None):
    """
        Function:
            spm_to_et
        Purpose:
            Calculates ephemeris time from a given SPM.
        Arguments:
            spm:
                Seconds-Past-Midnight.
            doy:
                The day of year.
            year:
                The given year as an integer.
        Keywords:
            kernels:
                Optional kernels to use for the geometry.
        Output:
            et:
                The ephemeris time from the given inputs.
    """

    if doy == 0:
        raise ValueError("Error: spm_to_et. input doy is 0!")

    if year == 0:
        raise ValueError("Error: spm_to_et. input year is 0!")

    # Leap seconds kernel
    if kernels is None:
        kernels = '../../kernels/naif/CASSINI/kernels/lsk/naif0012.tls'

    spiceypy.kclear()
    spiceypy.furnsh(kernels)

    hours = (spm / 3600.0).astype(int)
    minutes = ((spm - hours*3600.0).astype(int) / 60.0).astype(int)
    seconds = spm - hours*3600.0 - minutes*60.0
    n_pts = numpy.size(spm)
    et_sec_vals = numpy.zeros(n_pts)

    # Check for leap year
    if (year % 4) == 0:
        days_per_year = 366
    else:
        days_per_year = 365

    for ind in range(n_pts):

        # Adjust if hour makes it go on to next day
        doy_adjusted = doy + hours[ind]/24

        # Adjust if day goes on to next year
        this_year = (
            year + numpy.floor((doy_adjusted - 1) / days_per_year)
        ).astype(int)
        this_doy = (1 + (doy_adjusted - 1) % days_per_year).astype(int)

        utc_str = (
            str(this_year) + '-' + str(this_doy).zfill(3) + 'T' +
            str(hours[ind] % 24).zfill(2) + ':' +
            str(minutes[ind]).zfill(2) + ':' + f"{seconds[ind]:2.13F}"
        )

        et_sec_vals[ind] = spiceypy.utc2et(utc_str)

    return et_sec_vals
