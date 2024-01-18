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
#       Convert ephemeris time to seconds past midnight (SPM).                 #
################################################################################
#   Author: Jolene Fong                                                        #
#   Date:   2018/03/13                                                         #
################################################################################
"""

import pdb
import numpy
import spiceypy

def et_to_spm(et_vals, kernels = None, ref_doy = None):
    """
        Function:
            et_to_spm
        Purpose:
            Converts ephemeris time to seconds past midnight.
    """

    if kernels:
        spiceypy.kclear()
        spiceypy.furnsh(kernels)

    try:
        npts = len(et_vals)
    except TypeError:
        et_vals = [et_vals]
        npts = len(et_vals)

    spm_vals = []

    for ind in range(npts):

        # Convert ET to UTC string.
        utc_str = spiceypy.et2utc(et_vals[ind], "ISOD", 16)

        # Extract hour, minute, and second from UTC string.
        hour = float(utc_str[9:11])
        minute = float(utc_str[12:14])
        second = float(utc_str[15:])
        spm = hour*3600. + minute*60. + second
        spm_vals.append(spm)

    # Check if event goes through midnight.
    dr_vals = spm_vals - numpy.roll(spm_vals, 1)
    ibrk = numpy.argwhere(dr_vals[1:] < 0)

    if ibrk.size > 0:
        spm_vals_cont = numpy.asarray(spm_vals)

        for index in ibrk:
            ind = index[0]
            spm_brk = spm_vals[ind]
            spm_vals_cont[ind + 1:] = spm_vals_cont[ind+1:] + spm_brk

        spm_vals = spm_vals_cont

    if ref_doy:

        # Check if ref_doy is the same as current reference doy.
        utc_start = spiceypy.et2utc(et_vals[0], "ISOD", 16)
        doy_start = (utc_start[5:8])
        days_past = int(doy_start) - int(ref_doy)

        if days_past < 0:
            print("TROUBLE! Reference doy is after the first entry of et_vals!")
            pdb.set_trace()
        else:
            spm_vals = numpy.asarray(spm_vals) + days_past*24.*60.*60.

    return numpy.array(spm_vals)
