"""
    Script Name:
        rjm_rev7E_X43_e2e.py
    Purpose:
        Perform entire end-to-end process on rev7E X43 file taken
        from PDS node online.
    Revisions:
        2018/08/04 - rmaguire - Copied from gsteranka
        2018/08/06 - rmaguire - Edits to shell scripts.
"""

import numpy as np
import os
import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec
sys.path.append('..')
import rss_ringoccs as rss
sys.path.remove('..')
sys.path.append('../rss_ringoccs/tools')
from pds3_geo_series import write_geo_series
from pds3_cal_series import write_cal_series
from pds3_dlp_series import write_dlp_series
from pds3_tau_series import write_tau_series
sys.path.remove('../rss_ringoccs/tools')

rsr_file_name = 's10sroe2005123_0740nnnx43rd.2a2'
rsr_file_pds_dir = 'co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/'
rsr_file_local_dir = "../data/"+rsr_file_pds_dir
list_of_kernels = "../tables/rev007_list_of_kernels.txt"
kernels_list_file = '../tables/rev007_list_of_kernels.txt'
kernels_dir = '../kernels/'

print("Downloading RSR files...")
os.system('./get_rsr_file.sh %s %s %s ; echo "RSR Complete"' %
          (rsr_file_name, rsr_file_pds_dir, rsr_file_local_dir))
print("Downloading kernels...")
os.system('./get_kernels.sh %s %s ; echo "Kernels Complete"' %
          (kernels_list_file, kernels_dir))
rsr_file = rsr_file_local_dir + rsr_file_name
kernels = '../tables/Rev007_meta_kernel.ker'

output_directory = '../output/rev7E_X43_e2e_output/'
freq_offset_file = output_directory + 'freq_offset_file.txt'
f_resid_fit_parameters_file = output_directory + 'f_resid_fit_parameters.p'
power_norm_fit_parameters_file = (
    output_directory + 'power_norm_fit_parameters.p'
)
geo_file = 'RSS_2005_123_X43_E_GEO'
cal_file = 'RSS_2005_123_X43_E_CAL'
dlp_file = 'RSS_2005_123_X43_E_DLP'

f_USO = 8427222034.34050
dr_km_desired = 0.25
res_km = 1.0
inversion_range = [70000, 140000]
tau_file = 'RSS_2005_123_X43_E_TAU_' + str(int(res_km*1000)) + 'M'
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

write_geo_series(rev_info, geo_inst, geo_file, output_directory, 'Egress')
write_cal_series(rev_info, cal_inst, cal_file, output_directory, 'Egress')
write_dlp_series(rev_info, dlp_inst, dlp_file, output_directory, 'Egress')
write_tau_series(rev_info, tau_inst, tau_file, output_directory, 'Egress')


tau_data = -tau_inst.mu_vals*np.log(tau_inst.p_norm_vals)
# Make comparison plot.
plt.figure(figsize=(8.5, 11))
plt.suptitle("Rev007 E X43 - Maxwell Ringlet - 1000m Resolution",size=14)
gs = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0)

plt.subplot(gs[0,0])
plt.xlim(87400,87610)
plt.xticks((87410, 87460, 87510, 87560), ("87410", "87460", "87510", "87560"))
plt.ylim(-0.1,1.5)
plt.tick_params(
    axis='y', which='both', left=True, right=False,
    labelleft=True, labelright=False
)
plt.tick_params(
    axis='x', which='both', bottom=True, top=True,
    labelbottom=False, labeltop=False
)
plt.ylabel('Intensity')
plt.plot(tau_inst.rho_km_vals,tau_inst.p_norm_vals,label="Diffraction Profile")
plt.legend()

plt.subplot(gs[1,0])
plt.xlim(87400,87610)
plt.xticks((87410, 87460, 87510, 87560,87610), ("87410", "87460", "87510", "87560","87610"))
plt.ylim(-0.1,3.9)
plt.tick_params(
    axis='y', which='both', left=True, right=False,
    labelleft=True, labelright=False
)
plt.tick_params(
    axis='x', which='both', bottom=True, top=True,
    labelbottom=True, labeltop=False
)
plt.xlabel("Ring Radius (km)")
plt.ylabel('Optical Depth')
plt.plot(tau_inst.rho_km_vals,tau_data,label="Diffracted Optical Depth")
plt.legend()

plt.subplot(gs[0,1])
plt.xlim(87400,87610)
plt.xticks((87410, 87460, 87510, 87560), ("87410", "87460", "87510", "87560"))
plt.ylim(-0.1,1.5)
plt.tick_params(
    axis='y', which='both', left=False, right=True,
    labelleft=False, labelright=False
)
plt.tick_params(
    axis='x', which='both', bottom=True, top=True,
    labelbottom=False, labeltop=False
)
plt.plot(
    tau_inst.rho_km_vals,tau_inst.power_vals,
    label="Diffraction-Corrected Profile"
)
plt.legend()

plt.subplot(gs[1,1])
plt.xlim(87400,87610)
plt.xticks((87410, 87460, 87510, 87560), ("87410", "87460", "87510", "87560"))
plt.ylim(-0.1,3.9)
plt.tick_params(
    axis='y', which='both', left=False, right=True,
    labelleft=False, labelright=False
)
plt.tick_params(
    axis='x', which='both', bottom=True, top=True,
    labelbottom=True, labeltop=False
)
plt.xlabel("Ring Radius (km)")
plt.plot(tau_inst.rho_km_vals,tau_inst.tau_vals,label="Corrected Optical Depth")
plt.legend()
plt.show()
