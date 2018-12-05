
"""
pds3_reader.py

Purpose: Read legal PDS3 label files.

Revisions:
    2018 Jul 23 - jfong - copied from jwf_pds3_reader.py
    2018 Nov 05 - jfong - record all PDS keywords
"""
import pdb
import numpy as np

class PDS3Reader():
    def __init__(self, infile_lbl):

        # Read label file line by line
        line_list_lbl = []
        print('\nNow reading header:\n', infile_lbl)
        f = open(infile_lbl, 'r')
        for line in f:
            line_list_lbl.append(line)
        f.close()

        record_bytes, ind_rb = self.find_keyword_value(line_list_lbl,
                'RECORD_BYTES')
        file_records, ind_fr = self.find_keyword_value(line_list_lbl,
                'FILE_RECORDS')
        columns, ind_col = self.find_keyword_value(line_list_lbl, 'COLUMNS')
        frequency_band, ind_band = self.find_keyword_value(line_list_lbl, 
                'FREQUENCY_BAND')
        if frequency_band=='END':
            frequency_band, ind_band = self.find_keyword_value(line_list_lbl,
                    'BAND_NAME')


        obj_type, ind_obj = self.find_keyword_value(line_list_lbl, 'OBJECT')

        if obj_type != 'SERIES':
            print('TROUBLE! expect SERIES but found ', obj_type)
            pdb.set_trace()

        ncols = int(columns)
        names = ["" for x in range(ncols)]
        units = ["" for x in range(ncols)]
        contents = line_list_lbl
        line_found = ind_obj

        for n in range(ncols):
            contents = contents[line_found:]
            value, line_found = self.find_keyword_value(contents, 'OBJECT')
            if value != 'COLUMN':
                print('TROUBLE! expected COLUMN but found ', value)
                pdb.set_trace()

            contents    = contents[line_found:]
            name, line_found = self.find_keyword_value(contents, 'NAME')

            names[n] = name
            
            unit, line_found = self.find_keyword_value(contents, 'UNIT')

            units[n] = unit

            data_type, ind_data_type = self.find_keyword_value(contents,
                    'DATA_TYPE')
            if data_type != 'ASCII_REAL':
                print('TROUBLE! Cannot handle data_type = ', data_type)
                pdb.set_trace()
        series = SeriesReader(infile_lbl, names, record_bytes, file_records,
                frequency_band)

        self.series = series
        

    def find_keyword_value(self, line_list, keyword_to_find):
        nlines = len(line_list)
        for n in range(nlines):
            line = line_list[n]
            keyword, value = self.read_keyword_value(line)
            line_found = n + 1
            if keyword == keyword_to_find:
                break
        return value, line_found

    def read_keyword_value(self, line):
        eq_ind = line.find('=')

        keyword = line[:eq_ind].strip()
        value = line[eq_ind+1:].strip()

        return keyword, value

class SeriesReader():
    def __init__(self, infile_lbl, names, record_bytes, file_records, 
            frequency_band):

        infile_tab = str.replace(infile_lbl, 'LBL', 'TAB')
        product_id = infile_tab.split('/')[-1]

        self.RECORD_BYTES = record_bytes
        self.FILE_RECORDS = file_records
        self.PRODUCT_ID = product_id
        self.FREQUENCY_BAND = frequency_band

        column_list = []

        nlines = int(file_records)
        data_dict = {}
        # Initialize dictionary with lists... probs not the most efficient way
    
        for n in range(len(names)):
            data_dict[names[n]] = []

        print(infile_tab)
        print('\nNow reading table of data:\n', infile_tab)
        f = open(infile_tab, 'r')
        for line in f:
            columns = line.split(',')
            for col in range(len(columns)):
                val = float(columns[col].strip())
                data_dict[names[col]].append(val)

        f.close()

        self.convert_dict_to_attr(data_dict)
    
    def convert_dict_to_attr(self, input_dict):
        tc_naming_conventions = {
                'OBSERVED_EVENT_TIME': 't_oet_spm_vals'
                , 'RING_EVENT_TIME': 't_ret_spm_vals'
                , 'SPACECRAFT_EVENT_TIME': 't_set_spm_vals'
                , 'RING_RADIUS': 'rho_km_vals'
                , 'RING_LONGITUDE': 'phi_rl_deg_vals'
                , 'OBSERVED_RING_AZIMUTH': 'phi_ora_deg_vals'
                , 'OBSERVED_RING_ELEVATION': 'B_deg_vals'
                , 'SPACECRAFT_TO_RING_INTERCEPT_DISTANCE': 'D_km_vals'
                , 'RING_INTERCEPT_RADIAL_VELOCITY': 'rho_dot_kms_vals'
                , 'RING_INTERCEPT_AZIMUTHAL_VELOCITY': 'phi_rl_dot_kms_vals'
                , 'FRESNEL_SCALE': 'F_km_vals'
                , 'IMPACT_RADIUS': 'R_imp_km_vals'
                , 'SPACECRAFT_POSITION_X': 'rx_km_vals'
                , 'SPACECRAFT_POSITION_Y': 'ry_km_vals'
                , 'SPACECRAFT_POSITION_Z': 'rz_km_vals'
                , 'SPACECRAFT_VELOCITY_X': 'vx_kms_vals'
                , 'SPACECRAFT_VELOCITY_Y': 'vy_kms_vals'
                , 'SPACECRAFT_VELOCITY_Z': 'vz_kms_vals'
                , 'SKY_FREQUENCY': 'f_sky_hz_vals'
                , 'RESIDUAL_FREQUENCY': 'f_sky_resid_fit_vals'
                , 'FREE_SPACE_POWER': 'p_free_vals'
                , 'RADIUS_CORRECTION_DUE_TO_IMPROVED_POLE': 
                        'rho_corr_pole_km_vals'
                , 'RADIUS_CORRECTION_DUE_TO_TIMING_OFFSET':
                        'rho_corr_timing_km_vals'
                , 'NORMAL_OPTICAL_DEPTH': 'tau_vals'
                , 'NORMAL_OPTICAL_DEPTH_THRESHOLD': 'tau_threshold_vals'
                , 'PHASE_SHIFT': 'phase_deg_vals'
                }
        for key in input_dict:
            attr_name0 = str.replace(key.strip('"'), ' ', '_')
            attr_name = str.replace(attr_name0, '-', '_')
            attr_name_output = tc_naming_conventions[attr_name]
            attr_input = np.asarray(input_dict[key])

            setattr(self, attr_name_output, attr_input)
    

        return None

