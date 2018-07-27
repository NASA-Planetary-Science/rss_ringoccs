#!/usr/bin/env python
"""

calc_f_sky_recon.py

Purpose: Calculate sky frequency from the reconstructed event kernels.
         copied from Nicole Rappaport's "predicts" program in Fortran.

NOTE: Since this was copied directly with little idea of what is actually
      happening, the style looks horrible. Sorry about that

Revisions:
      gjs_calc_f_sky_recon.py
   2018 Feb 22 - gsteranka - Original version
      gjs_calc_f_sky_recon_v2.py
   2018 Mar 06 - gsteranka - Edited to be a set of functions instead of a
                             class, and to take SPM directly, instead of
                             ephemeris time
      calc_f_sky_recon.py
   2018 Mar 20 - gsteranka - Copy to official version and remove debug steps
"""

import numpy as np
from spiceypy import spiceypy as spice
import sys

try:
    from spm_to_et import spm_to_et
except ImportError:
    from ..tools.spm_to_et import spm_to_et

# Product of the gravitational constant by the Sun mass
MUS = 1.3271244002331E+11

# Parameter of the General Relativity Theory
GAMMA = 1.0000000000000E+00

# Factor (1+gamma)*mus/clight**3
LTfac = 9.8509818970128E-06

# see derpt.f for assignment of IDs
ID_SSB = [1, 2, 399, 4, 5, 6, 7, 8, 9, 10, 301, 606]

GM_SSB = [
    2.2032080000000E+04,  # 1
    3.2485859900000E+05,  #
    3.9860043600000E+05,  # 399
    4.2828314000000E+04,  # 4
    1.2671276786300E+08,  # 5
    3.7940626063000E+07,  # 6
    5.7945490070000E+06,  # 7
    6.8365340640000E+06,  # 8
    9.8160100000000E+02,  # 9
    1.3271244002331E+11,  # 10
    4.9027990000000E+03,  # 301
    8.9781370309840E+03]  # 606 (Titan) from cpck30Mar2016.tpc


def calc_f_sky_recon(f_spm, rsr_inst, sc_name, f_uso, kernels, TEST=False):
    """Primary function to calculate sky frequency at the requested times f_spm

    Args:
        f_spm (np.ndarray): SPM values to evaluate sky frequency at
        rsr_inst: Instance of RSRReader class
        sc_name (str): Name of spacecraft to get sky frequency for. In our
            case, this should always be 'Cassini'
        f_uso (float): USO sky frequency for the event and the right band
        kernels (list): String list of full path name to set of kernels
        TEST (bool): Print intermediate values if set to True
        """

    spice.kclear()
    spice.furnsh(kernels)

    et_vals = spm_to_et(f_spm, rsr_inst.doy, rsr_inst.year, kernels=kernels)

    # Spacecraft code
    sc_code = spice.bodn2c(sc_name)

    # Receiving Station code
    rs_code = spice.bodn2c(rsr_inst.dsn)

    n_times = len(et_vals)

    # Reconstructed sky frequency array
    RF = np.zeros(n_times)

    for i in range(n_times):
        et = et_vals[i]
        [etsc, _lt] = spice.ltime(et, rs_code, '<-', sc_code)
        A23 = derlt(sc_code, etsc, rs_code, et)
        B3 = derpt(et, rs_code)
        B2 = derpt(etsc, sc_code)
        temp = (B2 - B3)/(1.0 - B3)
        y = -A23*temp + A23 + temp
        y *= f_uso
        RF[i] = f_uso - y

    if TEST:
        for i in range(10):
            print('%30.16f' % RF[i])
    return RF


def derlt(sc_code, etsc, rs_code, et):
    """???"""

    ref = 'ECLIPJ2000'
    # Index for Solar System Barycenter
    SSB = 0
    abcorr = 'NONE'

    [S1, _lt] = spice.spkez(sc_code, etsc, ref, abcorr, SSB)
    [S2, _lt] = spice.spkez(rs_code, et, ref, abcorr, SSB)
    S12 = spice.vsubg(S2, S1, 6)

    R1 = spice.vnorm(S1[0:3])
    R2 = spice.vnorm(S2[0:3])
    R12 = spice.vnorm(S12[0:3])

    US1 = spice.vhat(S1[0:3])
    US2 = spice.vhat(S2[0:3])
    US12 = spice.vhat(S12[0:3])

    DR1DT1 = spice.vdot(US1, S1[3:6])
    DR2DT2 = spice.vdot(US2, S2[3:6])

    PR12T2 = spice.vdot(US12, S2[3:6])
    PR12T1 = -spice.vdot(US12, S1[3:6])

    D1 = R1 + R2 + R12
    D2 = R1 + R2 - R12

    TEMP1 = DR1DT1 + DR2DT2
    TEMP2 = PR12T1 + PR12T2

    TERM1 = (PR12T1 + PR12T2)/spice.clight()
    TERM2 = (TEMP1 + TEMP2)/D1 - (TEMP1 - TEMP2)/D2
    TERM2 *= LTfac
    NUMER = TERM2 + TERM1
    TERM1 = PR12T1/spice.clight()
    TERM2 = (DR1DT1 + PR12T1)/D1 - (DR1DT1 - PR12T1)/D2
    TERM2 *= LTfac
    DENOM = TERM2 + TERM1 + 1.0

    DLTDT2 = NUMER/DENOM
    return DLTDT2


def derpt(et, code):
    """???"""

    ref = 'ECLIPJ2000'
    # Index for Solar System Barycenter
    SSB = 0
    abcorr = 'NONE'

    [SI, _lt] = spice.spkez(code, et, ref, abcorr, SSB)
    SIDOT2 = np.sum(SI[3:6]**2)

    PHII = 0.0

    for i in range(len(ID_SSB)):
        body = ID_SSB[i]
        [SJ, _lt] = spice.spkez(body, et, ref, abcorr, SSB)
        SIJ = spice.vsubg(SJ, SI, 6)
        RIJ = spice.vnorm(SIJ[0:3])
        PHII += GM_SSB[i]/RIJ

    # From PRDCT_PARAM.PCK
    LSEC = 1.550520E-8
    B = (PHII + 0.5*SIDOT2)/(spice.clight()**2) - LSEC

    return B


if __name__ == '__main__':
    pass
