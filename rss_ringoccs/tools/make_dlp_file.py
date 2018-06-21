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
    2018 Jun 18 - gsteranka - Added spm_range keyword
"""

import numpy as np
import pdb
import spiceypy as spice

def make_dlp_file(norm_diff_inst, dlp_file_name, spm_range=None):
    """
    Write a DLP file from an instance of the associated class to a file of
    the specified name

    Args:
        norm_diff_inst: Instance of the NormDiff (normalized diffraction
            pattern) class
        dlp_file_name (str): File name to write DLP file to
        spm_range (list): minimum and maximum SPM value to save, in a list.
            If it's a chord occultation, make a list of 2 lists together,
            where the first list is ingress and the second list is egress:
            [[ing_min, ing_max], [egr_min, egr_max]]
    """

    fmt_str = ('%32.16f, %10i, %10i, ' + '%32.16f, '*4 + '%10i, '
        + '%32.16f, '*3 + '%32.16f')

    mu = np.sin(abs(norm_diff_inst.B_rad_vals))

    if norm_diff_inst.end_of_chord_ing is None:

        if spm_range is None:
            spm_range = [min(norm_diff_inst.spm_vals),
                max(norm_diff_inst.spm_vals)]
        ind = ((norm_diff_inst.spm_vals >= spm_range[0]) &
            (norm_diff_inst.spm_vals <= spm_range[1]))

        np.savetxt(dlp_file_name,
            np.c_[norm_diff_inst.rho_km_vals[ind],
                norm_diff_inst.rho_corr_pole_km_vals[ind],
                norm_diff_inst.rho_corr_timing_km_vals[ind],
                norm_diff_inst.phi_rl_rad_vals[ind]*spice.dpr(),
                norm_diff_inst.phi_rad_vals[ind]*spice.dpr(),
                -mu[ind]*np.log(norm_diff_inst.p_norm_vals[ind]),
                norm_diff_inst.phase_rad_vals[ind]*spice.dpr(),
                norm_diff_inst.tau_threshold_vals[ind],
                norm_diff_inst.spm_vals[ind],
                norm_diff_inst.t_ret_spm_vals[ind],
                norm_diff_inst.t_set_spm_vals[ind],
                norm_diff_inst.B_rad_vals[ind]*spice.dpr()],
            fmt=fmt_str)
    else:

        dlp_file_name_1 = dlp_file_name[0:-4] + '_I' + dlp_file_name[-4:]
        dlp_file_name_2 = dlp_file_name[0:-4] + '_E' + dlp_file_name[-4:]

        print('DETECTED CHORD OCCULTATION. SPLITTING INTO TWO FILES:')
        print(dlp_file_name_1 + '\n' + dlp_file_name_2)

        ind1 = norm_diff_inst.end_of_chord_ing
        ind2 = norm_diff_inst.end_of_chord_ing + 1

        if spm_range is None:
            spm_range = [[min(norm_diff_inst.spm_vals[0:ind1]),
                max(norm_diff_inst.spm_vals[0:ind1])],
                [min(norm_diff_inst.spm_vals[ind2:]),
                max(norm_diff_inst.spm_vals[ind2:])]]
        ind_ing = ((norm_diff_inst.spm_vals >= spm_range[0][0])
            & (norm_diff_inst.spm_vals <= spm_range[0][1]))
        ind_egr = ((norm_diff_inst.spm_vals >= spm_range[1][0])
            & (norm_diff_inst.spm_vals <= spm_range[1][1]))

        np.savetxt(dlp_file_name_1,
            np.c_[norm_diff_inst.rho_km_vals[ind_ing],
                norm_diff_inst.rho_corr_pole_km_vals[ind_ing],
                norm_diff_inst.rho_corr_timing_km_vals[ind_ing],
                norm_diff_inst.phi_rl_rad_vals[ind_ing]*spice.dpr(),
                norm_diff_inst.phi_rad_vals[ind_ing]*spice.dpr(),
                -mu[ind_ing]*np.log(norm_diff_inst.p_norm_vals[ind_ing]),
                norm_diff_inst.phase_rad_vals[ind_ing]*spice.dpr(),
                norm_diff_inst.tau_threshold_vals[ind_ing],
                norm_diff_inst.spm_vals[ind_ing],
                norm_diff_inst.t_ret_spm_vals[ind_ing],
                norm_diff_inst.t_set_spm_vals[ind_ing],
                norm_diff_inst.B_rad_vals[ind_ing]*spice.dpr()],
            fmt=fmt_str)

        np.savetxt(dlp_file_name_2,
            np.c_[norm_diff_inst.rho_km_vals[ind_egr],
                norm_diff_inst.rho_corr_pole_km_vals[ind_egr],
                norm_diff_inst.rho_corr_timing_km_vals[ind_egr],
                norm_diff_inst.phi_rl_rad_vals[ind_egr]*spice.dpr(),
                norm_diff_inst.phi_rad_vals[ind_egr]*spice.dpr(),
                -mu[ind_egr]*np.log(norm_diff_inst.p_norm_vals[ind_egr]),
                norm_diff_inst.phase_rad_vals[ind_egr]*spice.dpr(),
                norm_diff_inst.tau_threshold_vals[ind_egr],
                norm_diff_inst.spm_vals[ind_egr],
                norm_diff_inst.t_ret_spm_vals[ind_egr],
                norm_diff_inst.t_set_spm_vals[ind_egr],
                norm_diff_inst.B_rad_vals[ind_egr]*spice.dpr()],
            fmt=fmt_str)
