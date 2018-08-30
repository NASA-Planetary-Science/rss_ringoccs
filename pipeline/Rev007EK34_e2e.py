"""
Script Name:
    Rev007EK34_e2e.py
Purpose:
    Perform entire end-to-end process on rev7E K34 file taken
    from PDS node online. This is like the Rev007EX43_HuygensRinglet_test.py
    script, except there aren't files pre-saved for you, so you have to go
    through the steps yourself. The tutorial in
    ../tutorials/end_to_end_example.ipynb describes portions of this script
Notes:
    Frequency GUI inputs:
        
    Power GUI inputs:
        31057, 31064 ; 31783, 31790 ; 34074, 34082 ; 34240, 34252 ; 35275, 35295 ; 35540, 40000
        31061, 31787, 34078, 34246, 35285, 35550
        Fit Order = 1
Revisions:
    2018 Aug 04 - rmaguire - Copied from gsteranka
    2018 Aug 06 - rmaguire - Edits to shell scripts.
        gjs_rev7E_K34_e2e.py
    2018 Aug 21 - gsteranka - Separate version for rev7E K34, which doesn't
                              read or write pickle files. Made for the
                              purpose of GUI demonstration
        Rev007EK34_e2e.py
    2018 Aug 27 - gsteranka - Edited to include pickle file reading and writing
                              functions like in a normal full end-to-end run.
                              Also added GEO, CAL, DLP, and TAU file writing
                              functions.
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pdb
import pickle
import sys

sys.path.append('..')
import rss_ringoccs as rss
sys.path.remove('..')

sys.path.append('../rss_ringoccs/tools')
from pds3_geo_series import write_geo_series
from pds3_cal_series import write_cal_series
from pds3_dlp_series import write_dlp_series
from pds3_tau_series import write_tau_series
sys.path.remove('../rss_ringoccs/tools')

rsr_file_name = 's10sroe2005123_0740nnnk34rd.1b2'
rsr_file_pds_dir = 'co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/'
rsr_file_local_dir = "../data/" + rsr_file_pds_dir
kernels_list_file = '../tables/rev007_list_of_kernels.txt'
kernels_dir = '../kernels/'
kernels = '../tables/Rev007_meta_kernel.ker'
rev_number = '007'
occultation_direction = 'Egress'

print("Downloading RSR files...")
os.system('./get_rsr_file.sh %s %s %s ; echo "RSR Complete"' %
          (rsr_file_name, rsr_file_pds_dir, rsr_file_local_dir))
print("Downloading kernels...")
os.system('./get_kernels.sh %s %s ; echo "Kernels Complete"' %
          (kernels_list_file, kernels_dir))
rsr_file = rsr_file_local_dir + rsr_file_name

output_directory = '../output/rev7E_K34_e2e_output/'
freq_offset_file = output_directory + 'freq_offset_file.txt'
f_resid_fit_parameters_file = output_directory + 'f_resid_fit_parameters.p'
power_norm_fit_parameters_file = (
    output_directory + 'power_norm_fit_parameters.p'
)
geo_file = 'RSS_2005_123_K34_E_GEO'
cal_file = 'RSS_2005_123_K34_E_CAL'
dlp_file = 'RSS_2005_123_K34_E_DLP'

f_USO = 8427222034.34050 * 3.8
dr_km_desired = 0.25
res_km = 1.0
inversion_range = [87410, 87610]
tau_file = 'RSS_2005_123_K34_E_TAU'
verbose = True

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

rev_info = rss.tools.get_rev_info(rsr_inst, rev_number)

geo_inst = rss.occgeo.Geometry(
    rsr_inst, 'Saturn', 'Cassini', [kernels], verbose=verbose
)

os.system('[ ! -d ' + output_directory + ' ] && mkdir -p ' + output_directory)

# Calculate frequency offset if no file already there. Otherwise, read in the
#     previously made file
if os.path.exists(freq_offset_file):
    freq_offset_file_vals = np.loadtxt(freq_offset_file)
    f_spm = freq_offset_file_vals[:, 0]
    f_offset = freq_offset_file_vals[:, 1]
else:
    f_spm, f_offset, freq_offset_history = rss.calibration.calc_freq_offset(
        rsr_inst, freq_offset_file=freq_offset_file, verbose=verbose)

# Manually make fit to frequency offset if no file already there. Otherwise,
#     read in the fit parameters from the previously made file
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

# Manually make spline fit to unnormalized power if no file already there.
#     Otherwise, read in the fit parameters from the previously made file
if os.path.exists(power_norm_fit_parameters_file):
    k_power_norm, freespace_spm, knots_spm = read_power_norm_fit_parameters(
        power_norm_fit_parameters_file
    )
    spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
        freespace_spm=freespace_spm, knots_spm=knots_spm,
        spline_order=k_power_norm, USE_GUI=False, verbose=verbose
    )
else:
    spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(verbose=verbose)
    write_power_norm_fit_parameters(
        norm_inst, power_norm_fit_parameters_file
    )


cal_inst = rss.calibration.Calibration(
    fit_inst, norm_inst, geo_inst, verbose=verbose
)

dlp_inst = rss.calibration.NormDiff(
    rsr_inst, dr_km_desired, geo_inst, cal_inst, verbose=verbose
)

tau_inst = rss.diffcorr.DiffractionCorrection(
    dlp_inst, res_km, rng=inversion_range, verbose=verbose
)

write_geo_series(rev_info, geo_inst, geo_file, output_directory, occultation_direction)
write_cal_series(rev_info, cal_inst, cal_file, output_directory, occultation_direction)
write_dlp_series(rev_info, dlp_inst, dlp_file, output_directory, occultation_direction)
write_tau_series(rev_info, tau_inst, tau_file, output_directory, occultation_direction)

fig1 = plt.figure(1, figsize=(11, 8.5))
plt.plot(tau_inst.rho_km_vals, tau_inst.power_vals)
plt.title('Rev7E K34')
plt.xlabel('Rho (km)')
plt.ylabel('Normalized Power')
plt.show()
pdb.set_trace()
