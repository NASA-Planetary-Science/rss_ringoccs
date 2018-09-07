"""
    Script Name:
        gjs_power_normalization_tutorial_rev9E_X14_e2e.py
    Purpose:
        Demonstration of power normalization GUI
    Revisions:
        2018/08/04 - rmaguire - Copied from gsteranka
        2018/08/06 - rmaguire - Edits to shell scripts.
        2018/08/22/ - gsteranka - Copied to GUI demonstration
"""

import numpy as np
import os
import pdb
import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
sys.path.append('..')
import rss_ringoccs as rss
sys.path.remove('..')

rsr_file_name = 's11sroe2005159_1715nnnx14rd.2a2'
rsr_file_pds_dir = 'co-s-rss-1-sroc1-v10/cors_0112/sroc1_159/rsr/'
rsr_file_local_dir = '../data/' + rsr_file_pds_dir
kernels_list_file = '../tables/rev009_list_of_kernels.txt'
kernels_dir = '../kernels/'

print('Downloading RSR files...')
os.system('./get_rsr_file.sh %s %s %s ; echo "RSR Complete"' %
          (rsr_file_name, rsr_file_pds_dir, rsr_file_local_dir))
print('Downloading kernels...')
os.system('./get_kernels.sh %s %s ; echo "Kernels Complete"' %
          (kernels_list_file, kernels_dir))
rsr_file = rsr_file_local_dir + rsr_file_name
kernels = '../tables/Rev009_meta_kernel.ker'

output_directory = '../output/rev9E_X14_e2e_output/'
freq_offset_file = output_directory + 'freq_offset_file.txt'
f_resid_fit_parameters_file = output_directory + 'f_resid_fit_parameters.p'

f_USO = 8427222034.34050
dr_km_desired = 0.25
res_km = 1.0
inversion_range = [70000, 140000]
tau_file = 'RSS_2005_159_X14_E_TAU_' + str(int(res_km*1000)) + 'M'
verbose = True

# ***END OF USER INPUT**


def read_f_resid_fit_parameters(f_resid_fit_parameters_file):

    file_object = open(f_resid_fit_parameters_file, 'rb')
    fit_param_dict = pickle.load(file_object)
    k = fit_param_dict['k']
    spm_include = fit_param_dict['spm_include']
    return k, spm_include


def write_f_resid_fit_parameters(fit_inst, f_resid_fit_parameters_file):

    k = fit_inst._poly_order
    spm_include = fit_inst._spm_include
    fit_param_dict = {'k': k, 'spm_include': spm_include}
    file_object = open(f_resid_fit_parameters_file, 'wb')
    pickle.dump(fit_param_dict, file_object)
    file_object.close()


def read_power_norm_fit_parameters(power_norm_fit_parameters_file):

    file_object = open(power_norm_fit_parameters_file, 'rb')
    fit_param_dict = pickle.load(file_object)
    k = fit_param_dict['k']
    freespace_spm = fit_param_dict['freespace_spm']
    knots_spm = fit_param_dict['knots_spm']
    return k, freespace_spm, knots_spm


def write_power_norm_fit_parameters(norm_inst, power_norm_fit_parameters_file):

    k = norm_inst._spline_order
    freespace_spm = norm_inst._freespace_spm
    knots_spm = norm_inst._knots_spm
    fit_param_dict = {
        'k': k,
        'freespace_spm': freespace_spm,
        'knots_spm': knots_spm
    }
    file_object = open(power_norm_fit_parameters_file, 'wb')
    pickle.dump(fit_param_dict, file_object)
    file_object.close()


rsr_inst = rss.rsr_reader.RSRReader(rsr_file, verbose=verbose)

rev_info = rss.tools.get_rev_info(rsr_inst, '007')

geo_inst = rss.occgeo.Geometry(
    rsr_inst, 'Saturn', 'Cassini', [kernels], verbose=verbose
)

os.system('[ ! -d ' + output_directory + ' ] && mkdir -p ' + output_directory)

if os.path.exists(freq_offset_file):
    freq_offset_file_vals = np.loadtxt(freq_offset_file)
    f_spm = freq_offset_file_vals[:, 0]
    f_offset = freq_offset_file_vals[:, 1]
else:
    f_spm, f_offset, freq_offset_history = rss.calibration.calc_freq_offset(
        rsr_inst, freq_offset_file=freq_offset_file, verbose=verbose)

if os.path.exists(f_resid_fit_parameters_file):
    k_f_resid, spm_include = read_f_resid_fit_parameters(
        f_resid_fit_parameters_file
    )
    fit_inst = rss.calibration.FreqOffsetFit(
        rsr_inst, geo_inst, f_spm, f_offset, f_USO, poly_order=k_f_resid,
        spm_include=spm_include, USE_GUI=False, verbose=verbose
    )
else:
    fit_inst = rss.calibration.FreqOffsetFit(
        rsr_inst, geo_inst, f_spm, f_offset, f_USO, verbose=verbose
    )
    write_f_resid_fit_parameters(fit_inst, f_resid_fit_parameters_file)

spm_vals, IQ_c = fit_inst.get_IQ_c()

norm_inst = rss.calibration.Normalization(
    spm_vals, IQ_c, geo_inst, rsr_inst, verbose=verbose
)

spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(verbose=verbose)

cal_inst = rss.calibration.Calibration(
    fit_inst, norm_inst, geo_inst, verbose=verbose
)

dlp_inst = rss.calibration.NormDiff(
    rsr_inst, dr_km_desired, geo_inst, cal_inst, verbose=verbose
)

tau_inst = rss.diffcorr.DiffractionCorrection(
    dlp_inst, res_km, rng=inversion_range, verbose=verbose
)
