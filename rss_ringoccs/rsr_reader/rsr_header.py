"""

rsr_header.py

Purpose: Similar to RSRReader class in rsr_reader.py, but this only reads the
         header of an RSR file.

Revisions:
        rsr_header.py
    2018 Jun 11 - gsteranka - Original version
"""

import numpy as np
import os
import pdb
import struct
import sys


__endian = '>'

# Header contents of SFDU

__sfdu_field_names = [
    'sfdu_authority1', 'sfdu_authority2', 'sfdu_authority3',
    'sfdu_authority4',
    'sfdu_version',
    'sfdu_class',
    'sfdu_reserved1', 'sfdu_reserved2',
    'sfdu_data_desription1', 'sfdu_data_desription2',
    'sfdu_data_desription3', 'sfdu_data_desription4',
    'sfdu_length']
__sfdu_format = 'cccc'+'c'+'c'+'cc'+'cccc'+'Q'

# Header Aggregation
__ha_field_names = ['ha_type', 'ha_length']
__ha_format = 'H'+'H'

# Primary Header
__ph_field_names = [
    'ph_type',
    'ph_length',
    'ph_data_major',
    'ph_data_minor',
    'ph_mission_ID',
    'ph_format_code']
__ph_format = 'H'+'H'+'B'+'B'+'B'+'B'

# Secondary Header
__sh_field_names = [
    'sh_type',
    'sh_length',
    'sh_originator',
    'sh_last_modifier',
    'sh_rsr_software_id',
    'sh_record_sequence_number',
    'sh_spc_id',
    'sh_dss_id',
    'sh_rsr_id',
    'sh_schan_id',
    'sh_reserved',
    'sh_spacecraft',
    'sh_prdx_pass_number',
    'sh_ulband',
    'sh_dl_band',
    'sh_trk_mode',
    'sh_ul_dss_id',
    'sh_fgain_px_no',
    'sh_fgain_if_bandwidth',
    'sh_frov_flag',
    'sh_attenuation',
    'sh_adc_rms','sh_adc_peak',
    'sh_year','sh_doy','sh_seconds',
    'sh_bits_per_sample',
    'sh_data_error',
    'sh_sample_rate',
    'sh_ddc_lo',
    'sh_rfif_lo','sh_sfdu_year','sh_sfdu_doy','sh_sfdu_seconds',
    'sh_predicts_time_shift','sh_predicts_freq_override',
    'sh_predicts_freq_rate','sh_predicts_freq_offset',
    'sh_sub_channel_freq',
    'sh_rf_freq_point_1','sh_rf_freq_point_2','sh_rf_freq_point_3',
    'sh_schan_freq_point_1','sh_schan_freq_point_2',
    'sh_schan_freq_point_3',
    'sh_schan_freq_poly_coef_1','sh_schan_freq_poly_coef_2',
    'sh_schan_freq_poly_coef_3',
    'sh_schan_accum_phase',
    'sh_schan_phase_poly_coef_1','sh_schan_phase_poly_coef_2',
    'sh_schan_phase_poly_coef_3','sh_schan_phase_poly_coef_4',
    'sh_reserved2a','sh_reserved2b']
__sh_format = 'hh'+'BBh'+'hBB'+'BBcBHccBBbBBBBBHHIBBHHHHH'+22*'d'

# Data
__data_field_names = [
    'Data_type',
    'Data_length',
    'Data_QI']
__data_header_format = 'HH'

__field_names = (__sfdu_field_names + __ha_field_names + __ph_field_names
        + __sh_field_names + __data_field_names)

def rsr_header(rsr_file):
    """
    Purpose:
    Read header of RSR file, given just its full path name

    args:
        rsr_file (str): Full path name of RSR file
    """

    struct_hdr_fmt = (__endian + __sfdu_format + __ha_format
        + __ph_format + __sh_format
        + __data_header_format)
    struct_hdr_len = struct.calcsize(struct_hdr_fmt)
    struct_unpack_hdr = struct.Struct(struct_hdr_fmt).unpack_from

    try:
        with open(rsr_file, 'rb') as f:
            sfdu_hdr_raw = f.read(struct_hdr_len)
    except FileNotFoundError as err:
        print('ERROR (rsr_header): File not found! {}',format(err))
        sys.exit()

    # Unpack SFDU header
    sfdu_hdr = struct_unpack_hdr(sfdu_hdr_raw)
    sfdu_hdr_dict = dict(zip(__field_names, sfdu_hdr))

    # Find number of SFDU in file, and number of points per SFDU
    rsr_size = os.path.getsize(rsr_file)
    bytes_per_sfdu = sfdu_hdr_dict['sfdu_length'] + 20
    n_sfdu = rsr_size / bytes_per_sfdu
    sh_bits_per_sample = sfdu_hdr_dict['sh_bits_per_sample']
    bytes_per_sample = sh_bits_per_sample / 8
    data_length_per_sfdu = sfdu_hdr_dict['Data_length']
    n_pts_per_sfdu = np.int(data_length_per_sfdu / (2*bytes_per_sample))

    # Get array of SPM values for whole file
    sh_sample_rate_hz = sfdu_hdr_dict['sh_sample_rate'] * 1000.0
    sh_sfdu_seconds = sfdu_hdr_dict['sh_sfdu_seconds']
    dt = 1.0 / sh_sample_rate_hz
    end_spm_of_rsr = sh_sfdu_seconds + n_pts_per_sfdu*n_sfdu*dt
    n_pts = round((end_spm_of_rsr - sh_sfdu_seconds)/dt)
    spm_vals = float(sh_sfdu_seconds) + dt*np.arange(n_pts)

    out_dict = {'spm_vals': spm_vals, 'doy': sfdu_hdr_dict['sh_doy'],
        'year': sfdu_hdr_dict['sh_year'],
        'dsn': 'DSS-'+str(sfdu_hdr_dict['sh_dss_id']),
        'band': sfdu_hdr_dict['sh_dl_band'],
        'sample_rate_khz': sfdu_hdr_dict['sh_sample_rate']}

    if out_dict['year'] == 0:
        out_dict['year'] = sfdu_hdr_dict['sh_sfdu_year']
    if out_dict['doy'] == 0:
        out_dict['doy'] = sfdu_hdr_dict['sh_sfdu_doy']

    return out_dict


if __name__ == '__main__':
    rsr_file = '../../../data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'
    rsr_hdr = rsr_header(rsr_file)
    pdb.set_trace()
