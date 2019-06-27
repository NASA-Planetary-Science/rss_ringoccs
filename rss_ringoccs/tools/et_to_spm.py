"""

:Purpose: 
    Convert ephemeris time to seconds past midnight (SPM).

:Dependencies:
    #. numpy
    #. spiceypy
"""

import numpy as np
import spiceypy as spice

def et_to_spm(et_vals, kernels=None, ref_doy=None):
    """
    Convert ephemeris time to seconds past midnight.

    Arguments
        :et_vals (*float* or *np.ndarray*): ET seconds past J2000
    
    Keyword Arguments
        :kernels (*str* or *list*): Path to NAIF kernels
        :ref_doy (*int*): Reference day of year, typically used for
            occultations that occur over multiple days

    Returns
        :spm_vals (*float* or *np.ndarray*): Seconds past midnight
    """
    if kernels:
        spice.kclear()
        spice.furnsh(kernels)
    if isinstance(et_vals, float):
        npts = 1
        et_vals = [et_vals]
    else:
        npts = len(et_vals)
    if npts == 0:
        return []
        


    spm_vals = []

    for n in range(npts):
        et = et_vals[n]

        # Convert ET to UTC string
        utc_str = spice.et2utc(et, "ISOD", 16)
        
        # Extract hour, minute, and second from UTC string
        hour = float(utc_str[9:11])
        minute = float(utc_str[12:14])
        second = float(utc_str[15:])
        spm = hour*3600. + minute*60. + second
        spm_vals.append(spm)


    # Check if event goes through midnight
    dr_vals = spm_vals - np.roll(spm_vals, 1)
    ibrk = np.argwhere(dr_vals[1:] < 0)

    if ibrk.size > 0:
        spm_vals_cont = np.asarray(spm_vals)
        for index in range(len(ibrk)):
            ind = ibrk[index][0] 
            spm_brk = spm_vals[ind]
            spm_vals_cont[ind+1:] = spm_vals_cont[ind+1:] + spm_brk
        spm_vals = spm_vals_cont

    if ref_doy:
        # check if ref_doy is the same as current reference doy
        utc_start = spice.et2utc(et_vals[0], "ISOD", 16)
        doy_start = (utc_start[5:8])
        days_past = int(doy_start) - int(ref_doy)
        if days_past < 0:
            print("TROUBLE! Reference doy is after the first entry of et_vals!")
            pdb.set_trace()
        else:
            spm_vals = np.asarray(spm_vals) + days_past*24.*60.*60.
        

    return np.array(spm_vals)
