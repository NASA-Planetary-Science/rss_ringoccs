"""

RSS_loop_v2.py

Purpose: Initial example of script looping over RSR files with a meta-kernel.

Notes:
    [1] Meta-kernel won't work if it's not the same distance from kernels
        directory as this directory

Revisions:
        RSS_loop_v1.py
    2018 Jun 13 - gsteranka - Original version
        RSS_loop_v2.py
    2018 Jun 28 - gsteranka - Edited to be in new directory, and save to
                              directory structure matching Essam's
"""

import fnmatch
import numpy as np
import pandas as pd
import pdb
import pickle
import os
import sys

sys.path.append('..')
from rss_ringoccs.rsr_reader import RSRReader
from rss_ringoccs.rsr_reader import rsr_header
import rss_ringoccs.calibration as cal
import rss_ringoccs.occgeo as geo
import rss_ringoccs.tools as tools
#import rss_ringoccs.diffcorr as dc
sys.path.remove('..')

sys.path.append('../rss_ringoccs/diffcorr')
import gjs_diffraction_correction as dc
sys.path.remove('../rss_ringoccs/diffcorr')

sys.path.append('../rss_ringoccs/tools/')
from create_summary_doc import plot_summary_doc
from pds3_geo_series import write_geo_series
from pds3_cal_series import write_cal_series
from pds3_dlp_series import write_dlp_series
from pds3_tau_series import write_tau_series
from get_rev_info import get_rev_info
sys.path.remove('../rss_ringoccs/tools/')

planet = 'Saturn'
spacecraft = 'Cassini'

# CHANGE THIS FOR WHICH RSR FILES YOU WANT
# Change RSR files here
rsr_file_dir = '../../../../data/'
rsr_file_list = [#'s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1']#,
    #'s10-rev07-rsr-data/S10EAOE2005_123_0740NNNK34D.1B1']#,
    #'s10-rev07-rsr-data/S10EAOI2005_123_0230NNNX26D.3A1']#,
    #'s36-rev54-rsr-data/S36SROI2007353_0220NNNX63RD.1A1']#,
    #'s36-rev56-rsr-data/S36SROI2008015_2105NNNX63RD.1A1']#,
    #'cors_0112/sroc1_159/rsr/s11sroe2005159_1715nnnx14rd.2a2',
    #'cors_0113/sroc1_159/rsr/s11sroe2005159_1715nnnx55rd.1a2']
    'cors_0326/sroc12_169/rsr/S60SROE2010170_0330NNNK34RV.1N2']
rsr_files = [rsr_file_dir + rsr_file for rsr_file in rsr_file_list]

kernels = '../tables/Sa-TC17-V001.ker'

# X-band USO frequency
f_USO_X = 8427222034.34050

# Final spacing that you want in DLP and TAU files
dr_km_desired = 0.05

# CHANGE THIS FOR WHERE YOU WANT DIRECTORY TREES CREATED
# Directory in which to make new directories and output files
where_to_create_output = '../output/RSS_loop_v2_output/'

verbose = False


def date_to_rev(year, doy):
    """
    Pull rev number from a table given the year and doy. There are often
    multiple rows matching the year and doy, but these multiple rows always
    have the same rev number, so that's okay. Just take the first one.

    Args:
        year (int): Year of occultation
        doy (int): Day of year of occultation
    """

    date_to_rev_table = pd.read_csv(
        '../tables/RSSActivities_before_USOfailure_rings_only.txt',
        header=None, skiprows=1)
    rev_number_vals = np.asarray(date_to_rev_table[0])
    year_vals = np.asarray(date_to_rev_table[2])
    doy_vals = np.asarray(date_to_rev_table[3])
    year_where = (year_vals == year)
    doy_where = (doy_vals == doy)
    year_and_doy_where = year_where & doy_where
    rev_number = 'Rev' + rev_number_vals[year_and_doy_where][0][4:7]

    return rev_number


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
    fit_param_dict = {'k': k, 'freespace_spm': freespace_spm,
        'knots_spm': knots_spm}
    file_object = open(power_norm_fit_parameters_file, 'wb')
    pickle.dump(fit_param_dict, file_object)
    file_object.close()


# TODO (gsteranka): Test on both km and meters input
def res_km_str_to_float(res_km_str_vals):

    res_km_vals = []
    unit_vals = []
    for res_km_str in res_km_str_vals:
        if res_km_str[-2] == 'K':
            res_km_vals.append(float(res_km_str[:-2]))
            unit_vals.append('KM')
        else:
            res_km_vals.append(float(res_km_str[:-1])/1000.0)
            unit_vals.append('M')

    return res_km_vals, unit_vals


# Read Essam file information
EM_files_info_table = '../tables/EM_files_info_updated.txt'
EM_files_info = pd.read_csv(EM_files_info_table, delim_whitespace=True)
EM_files_info_rev_names = np.asarray(EM_files_info.RevName)
EM_files_info_tau_res = np.asarray(EM_files_info.tau_resolution)
EM_files_info_cal_spm_min = np.asarray(EM_files_info.cal_spm_min)
EM_files_info_cal_spm_max = np.asarray(EM_files_info.cal_spm_max)
EM_files_info_dlp_spm_min = np.asarray(EM_files_info.dlp_spm_min)
EM_files_info_dlp_spm_max = np.asarray(EM_files_info.dlp_spm_max)
EM_files_info_tau_spm_min = np.asarray(EM_files_info.tau_spm_min)
EM_files_info_tau_spm_max = np.asarray(EM_files_info.tau_spm_max)

# TODO (gsteranka): Uncomment this when ready to loop over all RSR files
# Get list of non-repeating 16kHz RSR files on our laptops. Ignore where it
#     says "no 16khz" files, and only take the first file where it gives a
#     2-element array of files, since they're always the same
# EM_files_info_rsr_16khz_list = EM_files_info.rsr_16khz_list.tolist()
# rsr_files = []
# for f in EM_files_info_rsr_16khz_list:
#     if ',' not in f:
#         f = f[2:-2]
#     elif 'no 16khz' in f:
#         #f = f[2:-2].split(':')[1]
#         continue
#     else:
#         f = f[2:-2].split(',')[0][:-1]
#     if ('../' + f) not in rsr_files:
#         rsr_files.append('../' + f)

print('Downloading kernels if needed')
os.system('cd ../kernels ; '
    + './get_kernels.sh ../tables/list_of_kernels.txt ; '
    +'cd ../pipeline')

for rsr_file in rsr_files:

    print('\nBeginning file: ' + rsr_file)

    # Read file to get year and DOY to determine what start/end SPM to
    #     use in EM_files_info
    rsr_hdr = rsr_header(rsr_file)
    rev_number = date_to_rev(rsr_hdr['year'], rsr_hdr['doy'])
    band = chr(rsr_hdr['band'][0])

    # Adjust USO frequency by wavelength
    if band == 'X':
        f_USO = f_USO_X
    elif band == 'S':
        f_USO = f_USO_X*(3.0/11.0)
    elif band == 'K':
        f_USO = f_USO_X*3.8
    else:
        sys.exit('ERROR: Band not implemented')

    print('Reading RSR file and decimating data')
    rsr_inst = RSRReader(rsr_file, decimate_16khz_to_1khz=True, verbose=verbose)
    rev_info = get_rev_info(rsr_inst, rev_number[3:6])

    print('Calculating geometry')
    geo_inst = geo.Geometry(rsr_inst, planet, spacecraft, [kernels])

    # Determine direction using rho_dot_kms_vals attribute of Geometry class
    if (np.any(geo_inst.rho_dot_kms_vals > 0)
            & np.any(geo_inst.rho_dot_kms_vals < 0)):
        direction = 'C'
        lbl_direction = ['Ingress', 'Egress']
    elif np.all(geo_inst.rho_dot_kms_vals > 0):
        direction = 'E'
        lbl_direction = 'Egress'
    elif np.all(geo_inst.rho_dot_kms_vals < 0):
        direction = 'I'
        lbl_direction = 'Ingress'
    else:
        sys.exit('ERROR: Bad GEO file. Check inputs to RSRReader and Geometry')

    rev_number_and_direction = rev_number + direction

    # Get start/end SPM of Essam's files
    rev_name_str = (rev_number_and_direction + '*_RSS_' + str(rsr_hdr['year'])
        + '_' + str(rsr_hdr['doy']).zfill(3) + '_' + band
        + rsr_hdr['dsn'][-2:] + '_*')
    rev_name_matches = fnmatch.filter(EM_files_info.RevName, rev_name_str)


    # SPM ranges for non-chord
    if direction != 'C':

        is_chord = False

        # If ingress/egress occultation, then all matches should be the same
        rev_name_match = rev_name_matches[0]

        row_vals = (EM_files_info_rev_names == rev_name_match)
        row_ind_vals = (np.arange(len(row_vals)))[row_vals]
        res_km_str_vals = EM_files_info_tau_res[row_vals].tolist()

        cal_spm_range_row = [EM_files_info_cal_spm_min[row_ind_vals[0]],
            EM_files_info_cal_spm_max[row_ind_vals[0]]]

        if not ((cal_spm_range_row[0] >= min(rsr_hdr['spm_vals'])) &
                (cal_spm_range_row[1] <= max(rsr_hdr['spm_vals']))):
            sys.exit('ERROR: Time of rev in table doesn\'t match actual times')

        print('Choosing spm range of Essam file ' + rev_name_match)
        dlp_spm_range_row = [EM_files_info_dlp_spm_min[row_ind_vals[0]],
            EM_files_info_dlp_spm_max[row_ind_vals[0]]]
        tau_spm_range_row = [EM_files_info.tau_spm_min[row_ind_vals[0]],
            EM_files_info.tau_spm_max[row_ind_vals[0]]]

    # SPM ranges for chord
    else:

        is_chord = True

        ingress_rev_name_matches = [s for s in rev_name_matches if 'I' in s]
        egress_rev_name_matches = [s for s in rev_name_matches if 'E' in s]

        rev_name_match_ing = ingress_rev_name_matches[0]
        rev_name_match_egr = egress_rev_name_matches[0]

        row_vals_ing = (EM_files_info_rev_names == rev_name_match_ing)
        row_vals_egr = (EM_files_info_rev_names == rev_name_match_egr)

        row_ind_vals_ing = (np.arange(len(row_vals_ing)))[row_vals_ing]
        row_ind_vals_egr = (np.arange(len(row_vals_egr)))[row_vals_egr]

        res_km_str_vals_ing = EM_files_info_tau_res[row_vals_ing].tolist()
        res_km_str_vals_egr = EM_files_info_tau_res[row_vals_egr].tolist()
        #res_km_str_vals = (res_km_str_vals_ing.tolist()
        #    + res_km_str_vals_egr.tolist())

        cal_spm_range_row = [EM_files_info_cal_spm_min[row_ind_vals_ing[0]],
            EM_files_info_cal_spm_max[row_ind_vals_ing[0]]]

        if not ((cal_spm_range_row[0] >= min(rsr_hdr['spm_vals'])) &
                (cal_spm_range_row[1] <= max(rsr_hdr['spm_vals']))):
            sys.exit('ERROR: Time of rev in table doesn\'t match actual times')

        print('Choosing spm range of Essam files ' + rev_name_match_ing
            + 'and ' + rev_name_match_egr)
        dlp_spm_range_ing = [EM_files_info_dlp_spm_min[row_ind_vals_ing[0]],
            EM_files_info_dlp_spm_max[row_ind_vals_ing[0]]]
        dlp_spm_range_egr = [EM_files_info_dlp_spm_min[row_ind_vals_egr[0]],
            EM_files_info_dlp_spm_max[row_ind_vals_egr[0]]]
        dlp_spm_range_row = [dlp_spm_range_ing, dlp_spm_range_egr]

        tau_spm_range_ing = [EM_files_info_tau_spm_min[row_ind_vals_ing[0]],
            EM_files_info_tau_spm_max[row_ind_vals_ing[0]]]
        tau_spm_range_egr = [EM_files_info_tau_spm_min[row_ind_vals_egr[0]],
            EM_files_info_tau_spm_max[row_ind_vals_egr[0]]]
        tau_spm_range_row = [tau_spm_range_ing, tau_spm_range_egr]


    # String pointing towards directory of current iteration's output files
    #     Also define files to look for or save for later
    if not is_chord:
        output_directory = (where_to_create_output + rev_number + '/'
            + direction + '/' + rev_number_and_direction + '_RSS_'
            + str(rsr_inst.year) + '_' + str(rsr_inst.doy).zfill(3) + '_'
            + band + rsr_inst.dsn[-2:] + '_' + direction + '/')
        freq_offset_file = output_directory + 'freq_offset_file.txt'
        f_resid_fit_parameters_file = (output_directory
            + 'f_resid_fit_parameters.p')
        power_norm_fit_parameters_file = (output_directory
            + 'power_norm_fit_parameters.p')
        geo_file = output_directory[-19:-1] + '_GEO'
        cal_file = output_directory[-19:-1] + '_CAL'
        dlp_file = output_directory[-19:-1] + '_DLP'

        # Convert resolution strings to integers
        res_km_vals, unit_vals = res_km_str_to_float(res_km_str_vals)

        tau_files = []
        for i in range(len(res_km_vals)):
            res_km = res_km_vals[i]
            _unit = unit_vals[i]
            if _unit == 'KM':
                tau_files.append(output_directory[-19:-1] + '_TAU_'
                    + str(int(res_km)).zfill(2) + _unit)
            elif _unit == 'M':
                tau_files.append(output_directory[-19:-1] + '_TAU_'
                    + str(int(res_km*1000)) + _unit)

        # Create output directory if it doesn't exist
        os.system('[ ! -d ' + output_directory + ' ] && mkdir -p '
            + output_directory)

    else:
        output_directory_I = (where_to_create_output + rev_number + '/'
            + 'I' + '/' + rev_number + 'CI' + '_RSS_' + str(rsr_inst.year)
            + '_' + str(rsr_inst.doy).zfill(3) + '_' + band + rsr_inst.dsn[-2:]
            + '_I/')
        output_directory_E = (where_to_create_output + rev_number + '/'
            + 'E' + '/' + rev_number + 'CE' + '_RSS_' + str(rsr_inst.year)
            + '_' + str(rsr_inst.doy).zfill(3) + '_' + band + rsr_inst.dsn[-2:]
            + '_E/')

        freq_offset_file_I = output_directory_I + 'freq_offset_file.txt'
        freq_offset_file_E = output_directory_E + 'freq_offset_file.txt'

        f_resid_fit_parameters_file_I = (output_directory_I
            + 'f_resid_fit_parameters.p')
        f_resid_fit_parameters_file_E = (output_directory_E
            + 'f_resid_fit_parameters.p')

        power_norm_fit_parameters_file_I = (output_directory_I
            + 'power_norm_fit_parameters.p')
        power_norm_fit_parameters_file_E = (output_directory_E
            + 'power_norm_fit_parameters.p')

        geo_file_I = output_directory_I[-19:-1] + '_GEO'
        geo_file_E = output_directory_E[-19:-1] + '_GEO'

        cal_file_I = output_directory_I[-19:-1] + '_CAL'
        cal_file_E = output_directory_E[-19:-1] + '_CAL'

        dlp_file_I = output_directory_I[-19:-1] + '_DLP'
        dlp_file_E = output_directory_E[-19:-1] + '_DLP'

        # Convert resolution strings to integers
        res_km_vals_ing, unit_vals_ing = res_km_str_to_float(res_km_str_vals_ing)
        res_km_vals_egr, unit_vals_egr = res_km_str_to_float(res_km_str_vals_egr)

        tau_files_I = []
        for i in range(len(res_km_vals_ing)):
            res_km = res_km_vals_ing[i]
            _unit = unit_vals_ing[i]
            if _unit == 'KM':
                tau_files_I.append(output_directory_I[-19:-1] + '_TAU_'
                    + str(int(res_km)).zfill(2) + _unit)
            elif _unit == 'M':
                tau_files_I.append(output_directory_I[-19:-1] + '_TAU_'
                    + str(int(res_km*1000)) + _unit)

        tau_files_E = []
        for i in range(len(res_km_vals_egr)):
            res_km = res_km_vals_egr[i]
            _unit = unit_vals_egr[i]
            if _unit == 'KM':
                tau_files_E.append(output_directory_E[-19:-1] + '_TAU_'
                    + str(int(res_km)).zfill(2) + _unit)
            elif _unit == 'M':
                tau_files_E.append(output_directory_E[-19:-1] + '_TAU_'
                    + str(int(res_km*1000)) + _unit)

        # Create output directories if they don't exist
        os.system('[ ! -d ' + output_directory_I + ' ] && mkdir -p '
            + output_directory_I)
        os.system('[ ! -d ' + output_directory_E + ' ] && mkdir -p '
            + output_directory_E)

    # Calculate frequency offset if there's no file to read from, otherwise
    #     read from the file already there
    print('Getting frequency offset')
    if not is_chord:
        if os.path.exists(freq_offset_file):
            freq_offset_file_vals = np.loadtxt(freq_offset_file)
            f_spm = freq_offset_file_vals[:, 0]
            f_offset = freq_offset_file_vals[:, 1]
        else:
            f_spm, f_offset, freq_offset_history = cal.calc_freq_offset(rsr_inst,
                freq_offset_file=freq_offset_file, verbose=verbose)
    else:
        if os.path.exists(freq_offset_file_I):
            freq_offset_file_vals = np.loadtxt(freq_offset_file_I)
            f_spm = freq_offset_file_vals[:, 0]
            f_offset = freq_offset_file_vals[:, 1]
        else:
            f_spm, f_offset, freq_offset_history = cal.calc_freq_offset(rsr_inst,
                freq_offset_file=freq_offset_file_I, verbose=verbose)
            os.system('cp ' + freq_offset_file_I + ' ' + freq_offset_file_E)

    # Use GUI to make fit to residual frequency if no file to read from,
    #     otherwise read and use parameters from file
    print('Getting residual frequency fit')
    if not is_chord:
        if os.path.exists(f_resid_fit_parameters_file):
            k_f_resid, spm_include = read_f_resid_fit_parameters(
                f_resid_fit_parameters_file)
            fit_inst = cal.FreqOffsetFit(rsr_inst, geo_inst, f_spm, f_offset,
                f_USO, poly_order=k_f_resid, spm_include=spm_include,
                USE_GUI=False, verbose=verbose)
        else:
            fit_inst = cal.FreqOffsetFit(rsr_inst, geo_inst, f_spm, f_offset,
                f_USO, verbose=verbose)
            write_f_resid_fit_parameters(fit_inst, f_resid_fit_parameters_file)
    else:
        if os.path.exists(f_resid_fit_parameters_file_I):
            k_f_resid, spm_include = read_f_resid_fit_parameters(
                f_resid_fit_parameters_file_I)
            fit_inst = cal.FreqOffsetFit(rsr_inst, geo_inst, f_spm, f_offset,
                f_USO, poly_order=k_f_resid, spm_include=spm_include,
                USE_GUI=False, verbose=verbose)
        else:
            fit_inst = cal.FreqOffsetFit(rsr_inst, geo_inst, f_spm, f_offset,
                f_USO, verbose=verbose)
            write_f_resid_fit_parameters(fit_inst, f_resid_fit_parameters_file_I)
            write_f_resid_fit_parameters(fit_inst, f_resid_fit_parameters_file_E)

    spm_vals, IQ_c = fit_inst.get_IQ_c()

    # Use GUI to make spline fit to power if no file to read from, otherwise
    #     read and use parameters from file
    print('Getting power normalization fit')
    norm_inst = cal.Normalization(spm_vals, IQ_c, geo_inst, rsr_inst,
        verbose=verbose)
    if not is_chord:
        if os.path.exists(power_norm_fit_parameters_file):
            k_power_norm, freespace_spm, knots_spm = read_power_norm_fit_parameters(
                power_norm_fit_parameters_file)
            spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
                freespace_spm=freespace_spm, knots_spm=knots_spm,
                spline_order=k_power_norm, USE_GUI=False, verbose=verbose)
        else:
            spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
                verbose=verbose)
            write_power_norm_fit_parameters(norm_inst,
                power_norm_fit_parameters_file)
    else:
        if os.path.exists(power_norm_fit_parameters_file_I):
            k_power_norm, freespace_spm, knots_spm = read_power_norm_fit_parameters(
                power_norm_fit_parameters_file_I)
            spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
                freespace_spm=freespace_spm, knots_spm=knots_spm,
                spline_order=k_power_norm, USE_GUI=False, verbose=verbose)
        else:
            spm_power_fit, power_spline_fit = norm_inst.get_spline_fit(
                verbose=verbose)
            write_power_norm_fit_parameters(norm_inst,
                power_norm_fit_parameters_file_I)
            write_power_norm_fit_parameters(norm_inst,
                power_norm_fit_parameters_file_E)

    cal_inst = cal.Calibration(fit_inst, norm_inst, geo_inst, verbose=verbose)

    print('Writing GEO TAB and LBL files')
    try:
        if not is_chord:
            write_geo_series(rev_info, geo_inst, geo_file, output_directory,
                lbl_direction)
        else:
            write_geo_series(rev_info, geo_inst, geo_file_I, output_directory_I,
                'Ingress')
            write_geo_series(rev_info, geo_inst, geo_file_E, output_directory_E,
                'Egress')
    except AttributeError:
        geo_inst.naif_toolkit_version = 'TEMPORARY'
        geo_inst.FREQUENCY_BAND = 'TEMPORARY'
        if not is_chord:
            write_geo_series(rev_info, geo_inst, geo_file, output_directory,
                lbl_direction)
        else:
            write_geo_series(rev_info, geo_inst, geo_file_I, output_directory_I,
                'Ingress')
            write_geo_series(rev_info, geo_inst, geo_file_E, output_directory_E,
                'Egress')

    print('Writing CAL TAB and LBL files')
    if not is_chord:
        write_cal_series(rev_info, cal_inst, cal_file, output_directory,
            lbl_direction)
    else:
        write_cal_series(rev_info, cal_inst, cal_file_I, output_directory_I,
            'Ingress')
        write_cal_series(rev_info, cal_inst, cal_file_E, output_directory_E,
            'Egress')

    print('Creating instance of NormDiff')
    norm_diff_inst = cal.NormDiff(rsr_inst, dr_km_desired, geo_inst, cal_inst,
        is_chord=is_chord, verbose=verbose)

    # TODO (gsteranka): Update chord section to use NormDiff.chord_split()
    print('Writing DLP file')
    if not is_chord:
        write_dlp_series(rev_info, norm_diff_inst, dlp_file, output_directory,
            lbl_direction)
    else:
        dlp_inst_I, dlp_inst_E = norm_diff_inst.chord_split()
        write_dlp_series(rev_info, dlp_inst_I, dlp_file_I, output_directory_I,
            'Ingress')
        write_dlp_series(rev_info, dlp_inst_E, dlp_file_E, output_directory_E,
            'Egress')

    # TODO (gsteranka): Include chord occultation and test it
    print('Performing diffraction correction')
    if not is_chord:
        for i in range(len(res_km_vals)):
            res_km = res_km_vals[i]
            _unit = unit_vals[i]
            tau_file = tau_files[i]
            dc_inst = dc.diffraction_correction(norm_diff_inst, res_km,
                rng=[70000, 150000])
            write_tau_series(rev_info, dc_inst, tau_file, output_directory,
                lbl_direction)
            if i == 0:
                print('Write summary file')
                plot_summary_doc(geo_inst, cal_inst, norm_diff_inst, dc_inst,
                    output_directory + rev_number_and_direction + '_' +
                    output_directory[-19:-1] + '_Summary.pdf')
    else:
        print('Chord diffraction correction not implemented yet!')

        for i in range(len(res_km_vals_ing)):
            res_km = res_km_vals_ing[i]
            _unit = unit_vals_ing[i]
            tau_file = tau_files_I[i]
            dc_inst = dc.diffraction_correction(dlp_inst_I, res_km,
                rng=[70000, 150000])
            write_tau_series(rev_info, dc_inst, tau_file, output_directory_I,
                'Ingress')
            if i == 0:
                print('Write summary file')
                plot_summary_doc(geo_inst, cal_inst, dlp_inst_I, dc_inst,
                    output_directory_I + rev_number + 'CI' + '_' +
                    output_directory_I[-19:-1] + '_Summary.pdf')

        for i in range(len(res_km_vals_egr)):
            res_km = res_km_vals_egr[i]
            _unit = unit_vals_egr[i]
            tau_file = tau_files_E[i]
            dc_inst = dc.diffraction_correction(dlp_inst_E, res_km,
                rng=[70000, 150000])
            write_tau_series(rev_info, dc_inst, tau_file, output_directory_E,
                'Egress')
            if i == 0:
                print('Write summary file')
                plot_summary_doc(geo_inst, cal_inst, dlp_inst_E, dc_inst,
                    output_directory_E + rev_number + 'CE' + '_' +
                    output_directory_E[-19:-1] + '_Summary.pdf')
