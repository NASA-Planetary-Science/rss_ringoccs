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
import numpy
import spiceypy
from .get_spm_from_et import get_spm_from_et

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

    npts = numpy.size(et_vals)
    spm_vals = [0] * npts

    # Compute the SPM values from the ET ones.
    for ind in range(npts):
        spm_vals[ind] = get_spm_from_et(et_vals[ind])

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
        doy_start = utc_start[5:8]
        days_past = int(doy_start) - int(ref_doy)

        if days_past < 0:
            raise ValueError("Reference DOY after first entry of et_vals.")

    return numpy.asarray(spm_vals) + days_past*86400.0
