"""

make_dlp_file.py

Purpose: Save a DLP file from an instance of the NormDiff class. These follow
         the same format as Essam

Revisions:
        make_dlp_file.py
    2018 May 23 - gsteranka - Original version
    2018 Jun 12 - gsteranka - Edited non-chord occultation section to work for
                              ingress. Didn't before b/c et_to_spm function
                              demands input times to be increasing. Originally
                              corrected this for chord occultations, but forgot
                              to do the same for non-chord
"""

import numpy as np
import pdb
import spiceypy as spice

try:
    from et_to_spm import et_to_spm
except ImportError:
    from .et_to_spm import et_to_spm


def make_dlp_file(norm_diff_inst, dlp_file_name):

    fmt_str = ('%32.16f, %10i, %10i, ' + '%32.16f, '*4 + '%10i, '
        + '%32.16f, '*3 + '%32.16f')

    mu = np.sin(abs(norm_diff_inst.B_rad_vals))

    # Filler arrays while we don't yet have these values
    n_pts = len(norm_diff_inst.rho_km_vals)
    TEMP = np.zeros(n_pts) - 999

    if norm_diff_inst.end_of_chord_ing is None:

        # et_to_spm function only takes spm increasing
        try:
            t_ret_spm_vals = et_to_spm(norm_diff_inst.t_ret_et_vals)
            t_set_spm_vals = et_to_spm(norm_diff_inst.t_set_et_vals)
        except TypeError:
            t_ret_spm_vals = (et_to_spm(norm_diff_inst.t_ret_et_vals[::-1]))[::-1]
            t_set_spm_vals = (et_to_spm(norm_diff_inst.t_set_et_vals[::-1]))[::-1]

        np.savetxt(dlp_file_name,
            np.c_[norm_diff_inst.rho_km_vals, TEMP, TEMP,
                norm_diff_inst.phi_rl_rad_vals*spice.dpr(),
                norm_diff_inst.phi_rad_vals*spice.dpr(),
                -mu*np.log(norm_diff_inst.p_norm_vals),
                norm_diff_inst.phase_rad_vals*spice.dpr(),
                TEMP, norm_diff_inst.spm_vals,
                t_ret_spm_vals, t_set_spm_vals,
                norm_diff_inst.B_rad_vals*spice.dpr()],
            fmt=fmt_str)
    else:

        dlp_file_name_1 = dlp_file_name[0:-4] + '_I' + dlp_file_name[-4:]
        dlp_file_name_2 = dlp_file_name[0:-4] + '_E' + dlp_file_name[-4:]

        print('DETECTED CHORD OCCULTATION. SPLITTING INTO TWO FILES:')
        print(dlp_file_name_1 + '\n' + dlp_file_name_2)

        ind1 = norm_diff_inst.end_of_chord_ing
        ind2 = norm_diff_inst.end_of_chord_ing + 1

        # et_to_spm function only takes spm increasing
        t_ret_et_vals_ing = norm_diff_inst.t_ret_et_vals[0:ind1]
        t_ret_et_vals_ing_reverse = t_ret_et_vals_ing[::-1]
        t_ret_spm_vals_ing_reverse = et_to_spm(t_ret_et_vals_ing_reverse)
        t_ret_spm_vals_ing = t_ret_spm_vals_ing_reverse[::-1]

        # et_to_spm function only takes spm increasing
        t_set_et_vals_ing = norm_diff_inst.t_set_et_vals[0:ind1]
        t_set_et_vals_ing_reverse = t_set_et_vals_ing[::-1]
        t_set_spm_vals_ing_reverse = et_to_spm(t_set_et_vals_ing_reverse)
        t_set_spm_vals_ing = t_set_spm_vals_ing_reverse[::-1]

        np.savetxt(dlp_file_name_1,
            np.c_[norm_diff_inst.rho_km_vals[0:ind1], TEMP[0:ind1],
                TEMP[0:ind1],
                norm_diff_inst.phi_rl_rad_vals[0:ind1]*spice.dpr(),
                norm_diff_inst.phi_rad_vals[0:ind1]*spice.dpr(),
                -mu[0:ind1]*np.log(norm_diff_inst.p_norm_vals[0:ind1]),
                norm_diff_inst.phase_rad_vals[0:ind1]*spice.dpr(),
                TEMP[0:ind1], norm_diff_inst.spm_vals[0:ind1],
                t_ret_spm_vals_ing, t_set_spm_vals_ing,
                norm_diff_inst.B_rad_vals[0:ind1]*spice.dpr()],
            fmt=fmt_str)

        np.savetxt(dlp_file_name_2,
            np.c_[norm_diff_inst.rho_km_vals[ind2:], TEMP[ind2:], TEMP[ind2:],
                norm_diff_inst.phi_rl_rad_vals[ind2:]*spice.dpr(),
                norm_diff_inst.phi_rad_vals[ind2:]*spice.dpr(),
                -mu[ind2:]*np.log(norm_diff_inst.p_norm_vals[ind2:]),
                norm_diff_inst.phase_rad_vals[ind2:]*spice.dpr(),
                TEMP[ind2:], norm_diff_inst.spm_vals[ind2:],
                et_to_spm(norm_diff_inst.t_ret_et_vals[ind2:]),
                et_to_spm(norm_diff_inst.t_set_et_vals[ind2:]),
                norm_diff_inst.B_rad_vals[ind2:]*spice.dpr()],
            fmt=fmt_str)
