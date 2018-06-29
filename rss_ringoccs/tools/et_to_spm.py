
"""

et_to_spm.py

Purpose: Convert ephemeris time to seconds past midnight (SPM).

Revisions:
    2018 Mar 13 - jfong - original
    2018 Jun 26 - jfong - add load_kernels if kernels present
                        - debug for occs longer than a day
    2018 Jun 28 - jfong - add ref_doy keyword; if set, calculate spm relative
                          to ref_doy
"""

import numpy as np
import spiceypy as spice
import pdb

def et_to_spm(et_vals, kernels=None, ref_doy=None):
    """Function to convert ephemeris time to SPM"""

    if kernels:
        spice.kclear()
        for kernel in kernels:
            spice.furnsh(kernel)

    npts = len(et_vals)

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
        print(utc_str, hour, minute, second, spm)
        spm_vals.append(spm)

    # Check if event goes through midnight
    dr_vals = spm_vals - np.roll(spm_vals, 1)
#    ibrk = np.argwhere(dr_vals[1:] < 0)
    ibrk = np.argwhere(dr_vals[1:] < 0)

    if ibrk.size > 0:
        spm_vals_cont = np.asarray(spm_vals)
        for index in range(len(ibrk)):
            ind = ibrk[index][0] 
            spm_brk = spm_vals[ind-1]
            spm_vals_cont[ind:] = spm_vals_cont[ind:] + spm_brk
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
            spm_vals = spm_vals + days_past*24.*60.*60.
        

    return np.array(spm_vals)


