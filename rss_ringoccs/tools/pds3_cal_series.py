
"""
pds3_cal_series.py

Purpose: From a cal instance, produce inputs to pds3_write_series_lbl()
         for *CAL.LBL

Revisions:
    2017 Jul 23 - jfong - copied from jwf_pds3_cal_series_v2.py
"""

import pds3_write_series as pds3
import numpy as np
import pdb
import time

def write_cal_series_data(cal_inst, fmt, out_file):
    """
    This writes a CAL data file.

    Args:
        cal_inst (class): Instance of Calibration class
        fmt (str): Format string
        out_file (str): Output file name, including path.
    """
    fmt_comma = fmt+','
    format_str = fmt_comma*3 + fmt + '%s'
    npts = len(cal_inst.t_oet_spm_vals)

    print('\nWriting CAL data to: ', out_file, '\n')
    f = open(out_file, 'w')
    for n in range(npts):
        f.write(format_str % (
            cal_inst.t_oet_spm_vals[n],
            cal_inst.f_sky_hz_vals[n],
            cal_inst.f_sky_resid_fit_vals[n],
            cal_inst.p_free_vals[n],
            '\r\n'))

    f.close()

    return None


def get_cal_series_info(rev_info, cal_inst, series_name, prof_dir):
    """
    This returns the information needed to write a CAL label file.

    Args:
        rev_info (dict): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        cal_inst (class): Instance of Calibration class
        series_name (str): Name of the output .TAB and .LBL file, not including
                           extensions. Date in YYYYMMDD format will be added
                           onto series_name
        prof_dir (str): Direction of ring occultation for this cal_inst

    Outputs:
        str_lbl (dict): Dictionary with keys: string_delimiter,
                        alignment_column, series_alignment_column,
                        keywords_value, keywords_NAIF_TOOLKIT_VERSION,
                        description, keywords_series, object_keys,
                        object_values, history

    Notes:
        [1] This is a reproduction of CAL label files within
            Cassini_RSS_Ring_Profiles_2018_Archive, with minor edits.
        [2] The format of each data entry is hardcoded within "nchar".
    """
    # Get current time in _YYYYMMDD-hhmmss format
    current_time_ISOD = time.strftime("%Y-%j") + 'T' + time.strftime("%H:%M:%S")

    # Values for number of columns, number of bytes per column,
    #   number of bytes per column delimiter, number of bytes allocated to
    #   special characters per column
    ncol = 4        
    nchar = 32      
    ndelim = 1      
    nspecial = 2    

    # Values for aligning equal signs
    alignment_column        = 32
    series_alignment_column = 28
    bytes_list = [nchar] * ncol

    # Calculate start bytes for each column
    new_bytes_list = [(nchar+ndelim)] * ncol
    new_bytes_list.insert(0,1)
    start_bytes_list = np.cumsum(new_bytes_list)
    start_bytes_list = start_bytes_list[:-1]

    col_num = list(range(1, ncol+1))

    rsr_file = rev_info['rsr_file']
    band = rev_info['band']
    year = rev_info['year']
    doy = rev_info['doy']
    dsn = rev_info['dsn']

    # Extract relevant information from cal instance
    sampling_parameter_arr = cal_inst.t_oet_spm_vals
    t_oet_spm_start = cal_inst.t_oet_spm_vals[0]
    t_oet_spm_end = cal_inst.t_oet_spm_vals[-1]
    geo_kernels = cal_inst.history['Input Variables']['geo_inst'][
            'Input Variables']['kernels']
    # Remove directory path from kernel list
    geo_kernels = ['"'+x.split('/')[-1]+'"' for x in geo_kernels]


    PDS_VERSION_ID = 'PDS3'
    RECORD_TYPE = 'FIXED LENGTH'
    RECORD_BYTES = pds3.get_record_bytes(ncol, nchar, ndelim, nspecial) 
    FILE_RECORDS = str(len(sampling_parameter_arr))
    SERIES_NAME = series_name

    DATA_SET_ID = '"CO-SR-RSS-?/?-OCC-V0.1"'
    RING_OBSERVATION_ID = pds3.get_ring_obs_id(year, doy, band, dsn)
    PRODUCT_ID = series_name
    PRODUCT_TYPE = 'CALIBRATION PARAMETERS'
    PRODUCT_CREATION_TIME = current_time_ISOD
    PRODUCER_ID = '"TC2017"'
    SOURCE_PRODUCT_ID = '"' + rsr_file.upper() + '"'


    INSTRUMENT_HOST_NAME = '"CASSINI ORBITER"'
    INSTRUMENT_HOST_ID = 'CO'
    INSTRUMENT_NAME = '"RADIO SCIENCE SUBSYSTEM"'
    INSTRUMENT_ID = 'RSS'
    MISSION_PHASE_NAME = '"TOUR"'
    TARGET_NAME = '"S RINGS"'
    START_TIME = pds3.get_ISOD_str(t_oet_spm_start, year, doy)
    STOP_TIME = pds3.get_ISOD_str(t_oet_spm_end, year, doy)
    REVOLUTION_NUMBER = rev_info['rev_num']
    DSN_STATION_NUMBER = dsn.split('-')[-1]
    OCCULTATION_TYPE = 'RADIO'
    PLANETARY_OCCULTATION_FLAG = rev_info['planetary_occ_flag']
    RING_OCCULTATION_DIRECTION = rev_info['occ_dir']
    RING_PROFILE_DIRECTION = prof_dir 
    FREQUENCY_BAND = band

    # '' indicates that this keyword is not present in LBL file
    NAIF_TOOLKIT_VERSION = ''

    SPICE_FILE_NAME = geo_kernels






    PDS_VERSION_ID_dict = {
            'key_order': ['PDS_VERSION_ID', 'RECORD_TYPE',
                        'RECORD_BYTES', 'FILE_RECORDS', '^SERIES']
            , 'PDS_VERSION_ID': PDS_VERSION_ID
            , 'RECORD_TYPE': RECORD_TYPE
            , 'RECORD_BYTES': RECORD_BYTES
            , 'FILE_RECORDS': FILE_RECORDS
            , '^SERIES': SERIES_NAME
            }

    DATA_SET_ID_dict = {
            'key_order': ['DATA_SET_ID', 'RING_OBSERVATION_ID',
                        'PRODUCT_ID', 'PRODUCT_TYPE',
                        'PRODUCT_CREATION_TIME', 'PRODUCER_ID',
                        'SOURCE_PRODUCT_ID']
            , 'DATA_SET_ID': DATA_SET_ID
            , 'RING_OBSERVATION_ID': RING_OBSERVATION_ID
            , 'PRODUCT_ID': PRODUCT_ID
            , 'PRODUCT_TYPE': PRODUCT_TYPE
            , 'PRODUCT_CREATION_TIME': PRODUCT_CREATION_TIME
            , 'PRODUCER_ID': PRODUCER_ID
            , 'SOURCE_PRODUCT_ID': SOURCE_PRODUCT_ID
            }

    INSTRUMENT_HOST_NAME_dict = {
            'key_order': ['INSTRUMENT_HOST_NAME'
                         , 'INSTRUMENT_HOST_ID'
                         , 'INSTRUMENT_NAME'
                         , 'INSTRUMENT_ID'
                         , 'MISSION_PHASE_NAME'
                         , 'TARGET_NAME'
                         , 'START_TIME'
                         , 'STOP_TIME'
                         , 'REVOLUTION_NUMBER'
                         , 'DSN_STATION_NUMBER'
                         , 'OCCULTATION_TYPE'
                         , 'PLANETARY_OCCULTATION_FLAG'
                         , 'RING_OCCULTATION_DIRECTION'
                         , 'RING_PROFILE_DIRECTION'
                         , 'FREQUENCY_BAND']
            , 'INSTRUMENT_HOST_NAME': INSTRUMENT_HOST_NAME
            , 'INSTRUMENT_HOST_ID': INSTRUMENT_HOST_ID
            , 'INSTRUMENT_NAME' : INSTRUMENT_NAME
            , 'INSTRUMENT_ID': INSTRUMENT_ID
            , 'MISSION_PHASE_NAME': MISSION_PHASE_NAME
            , 'TARGET_NAME': TARGET_NAME
            , 'START_TIME' : START_TIME
            , 'STOP_TIME' : STOP_TIME
            , 'REVOLUTION_NUMBER': REVOLUTION_NUMBER
            , 'DSN_STATION_NUMBER': DSN_STATION_NUMBER
            , 'OCCULTATION_TYPE' : OCCULTATION_TYPE
            , 'PLANETARY_OCCULTATION_FLAG': PLANETARY_OCCULTATION_FLAG
            , 'RING_OCCULTATION_DIRECTION': RING_OCCULTATION_DIRECTION
            , 'RING_PROFILE_DIRECTION': RING_PROFILE_DIRECTION
            , 'FREQUENCY_BAND': FREQUENCY_BAND
            }
    
    NAIF_TOOLKIT_VERSION_dict = {
            'key_order': ['NAIF_TOOLKIT_VERSION','SPICE_FILE_NAME']
            , 'NAIF_TOOLKIT_VERSION': NAIF_TOOLKIT_VERSION
            , 'SPICE_FILE_NAME': SPICE_FILE_NAME
            }


    qq = "'"

    sd = '|'
    FILE_DESCRIPTION = ('"This file contains estimates of' + sd
            + 'signal attributes needed to calibrate the raw ring data '
            + 'before' + sd + 'reliable optical depth profiles are '
            + 'computed. The attributes are the' + sd + 'signal sky-frequency, '
            + 'the fit to residual frequency, and the free-space' + sd
            + 'signal power. The attributes are listed versus ' 
            + 'OBSERVED EVENT TIME' + sd + '(Earth received time) ' 
            + 'over equal time increments (1 s).' + sd
            + ' ' + sd
            + 'The sky frequency estimates are included for completeness ' 
            + 'and to facilitate' + sd + ' independent checks of frequency '
            + 'calculations. The radio science receiver' + sd + 'at the DSN '
            + 'ground receiving station (the RSR) steers the frequency of the'
            + sd + 'received sinusoid so that the measured spectral '
            + 'line falls at the center of' + sd + 'the recording '
            + 'bandwidth. Spectral estimates of the measured I/Q samples can '
            + sd + 'be used to calculate any offsets of the spectral line '
            + 'from the center of' + sd + 'the bandwidth. The measured offsets '
            + 'together with other frequency' + sd + 'steering information '
            + 'encoded in the RSR recording are used to calculate the' + sd
            + 'listed sky-frequency based on procedures documented '
            + 'in JPLD-16765.' + sd
            + ' ' + sd
            + 'The frequency residual estimates are required to steer the '
            + 'frequency of the' + sd + 'downlink sinusoid to a constant '
            + 'value (here the center of the recording' + sd
            + 'bandwidth) before the sinusoid amplitude and phase can '
            + 'be estimated. This is' + sd + 'done in the FreqOffsetFit class '
            + '(inputs to this class are under "fit_inst"' + sd
            + 'history below.' + sd
            + ' ' + sd
            + 'The free-space signal power estimates are required to '
            + 'normalize the' + sd + 'power of the steered sinusoid to a '
            + 'value of about unity (optical depth' + sd + 'of about zero) '
            + 'outside Ring A, inside Ring C, and within large ring' + sd
            + 'gaps. This is done in the Normalization class (inputs are '
            + 'under "norm_inst"' + sd + 'history below).' + sd
            + ' ' + sd
            + 'This file was produced using the rss_ringoccs open-source '
            + 'processing suite' + sd + 'developed at Wellesley College with '
            + 'the support of the Cassini project and' + sd + 'hosted '
            + 'on GithHub at https://github.com/NASA-Planetary-Science/'
            + 'rss_ringoccs.' + sd
            + ' ' + sd
            + 'Please address any inquiries to:' + sd
            + 'Richard G. French' + sd
            + 'Astronomy Department, Wellesley College' + sd
            + 'Wellesley, MA 02481-8203' + sd
            + '(781) 283-3747' + sd
            + 'rfrench@wellesley.edu"')



    HIST_USER_NAME = cal_inst.history['User Name']
    HIST_HOST_NAME = cal_inst.history['Host Name']
    HIST_RUN_DATE = cal_inst.history['Run Date']
    HIST_PYTHON_VERSION = cal_inst.history['Python Version']
    HIST_OPERATING_SYSTEM = cal_inst.history['Operating System']
    HIST_SOURCE_DIR = cal_inst.history['Source Directory']
    HIST_SOURCE_FILE = cal_inst.history['Source File']
    HIST_INPUT_VARIABLES = cal_inst.history['Input Variables']
    HIST_INPUT_KEYWORDS = cal_inst.history['Input Keywords']
    HIST_description = ('This is a detailed record of the' + sd
                    + 'processing steps used to generate this file.')

    HISTORY_dict = {
            'key_order': ['User Name', 'Host Name', 'Run Date',
                        'Python Version', 'Operating System', 
                        'Source Directory','Source File', 
                        'Input Variables', 'Input Keywords']
            , 'hist name': 'Calibration history'
            , 'User Name': HIST_USER_NAME
            , 'Host Name': HIST_HOST_NAME
            , 'Run Date': HIST_RUN_DATE
            , 'Python Version': HIST_PYTHON_VERSION
            , 'Operating System': HIST_OPERATING_SYSTEM
            , 'Source Directory': HIST_SOURCE_DIR
            , 'Source File': HIST_SOURCE_FILE
            , 'Input Variables': HIST_INPUT_VARIABLES
            , 'Input Keywords': HIST_INPUT_KEYWORDS
            , 'description': HIST_description
            }


    blank = 'BLANK'
    keywords_values = [
        PDS_VERSION_ID_dict,
        blank,
        DATA_SET_ID_dict,
        blank,
        INSTRUMENT_HOST_NAME_dict,
        blank
        ]

    SERIES = 'SERIES'
    SERIES_NAME = '"CALIBRATION PARAMETERS"'
    SERIES_INTERCHANGE_FORMAT = 'ASCII'
    SERIES_COLUMNS = str(ncol)
    SERIES_ROWS = FILE_RECORDS
    SERIES_ROW_BYTES = RECORD_BYTES
    SERIES_SAMPLING_PARAMETER_NAME = '"OBSERVED EVENT TIME"'
    SERIES_SAMPLING_PARAMETER_UNIT = '"SECOND"'
    SERIES_MINIMUM_SAMPLING_PARAMETER = str(min(sampling_parameter_arr))
    SERIES_MAXIMUM_SAMPLING_PARAMETER = str(max(sampling_parameter_arr))
    SERIES_SAMPLING_PARAMETER_INTERVAL = pds3.get_sampling_interval(
            sampling_parameter_arr)
    SERIES_DESCRIPTION = ('"This series contains variations of the' + sd
            + 'signal sky-frequency, residual frequency, and free-space '
            + 'power as a' + sd + 'function of OBSERVED EVENT TIME '
            + '(Earth receiving time)."')
    
    OBJECT_REFERENCE_TIME = pds3.get_spm_ref_time(year, doy)


    SERIES_dict = {
            'key_order': ['OBJECT', 'NAME', 'INTERCHANGE_FORMAT',
                        'COLUMNS', 'ROWS', 'ROW_BYTES',
                        'SAMPLING_PARAMETER_NAME',
                        'SAMPLING_PARAMETER_UNIT',
                        'MINIMUM_SAMPLING_PARAMETER',
                        'MAXIMUM_SAMPLING_PARAMETER',
                        'SAMPLING_PARAMETER_INTERVAL',
                        'DESCRIPTION']
            ,'OBJECT': SERIES
            , 'NAME': SERIES_NAME
            , 'INTERCHANGE_FORMAT': SERIES_INTERCHANGE_FORMAT
            , 'COLUMNS': SERIES_COLUMNS
            , 'ROWS': SERIES_ROWS
            , 'ROW_BYTES': SERIES_ROW_BYTES
            , 'SAMPLING_PARAMETER_NAME': SERIES_SAMPLING_PARAMETER_NAME
            , 'SAMPLING_PARAMETER_UNIT': SERIES_SAMPLING_PARAMETER_UNIT
            , 'MINIMUM_SAMPLING_PARAMETER':
                SERIES_MINIMUM_SAMPLING_PARAMETER
            , 'MAXIMUM_SAMPLING_PARAMETER':
                SERIES_MAXIMUM_SAMPLING_PARAMETER
            , 'SAMPLING_PARAMETER_INTERVAL':
                SERIES_SAMPLING_PARAMETER_INTERVAL
            , 'DESCRIPTION': SERIES_DESCRIPTION
            }

    object_keys = [
            'NAME'
            , 'COLUMN_NUMBER'
            , 'DATA_TYPE'
            , 'START_BYTE'
            , 'BYTES'
            , 'FORMAT'
            , 'UNIT'
            , 'REFERENCE_TIME'
            , 'DESCRIPTION'
            ]

    object_names = [
            '"OBSERVED EVENT TIME"'
            , '"SKY FREQUENCY"'
            , '"RESIDUAL FREQUENCY"'
            , '"FREESPACE POWER"'
            ]

    n_objects = len(object_names)
    data_types = ['ASCII_REAL'] * n_objects
    formats = ['"F32.16"'] * n_objects
    units = ['"SECOND"', '"HERTZ"', '"HERTZ"', '"N/A"']
    
    es = ''
    reference_times = [OBJECT_REFERENCE_TIME, es, es, es]


    object_descriptions = [
            ('"The instant at which photons were' + sd 
            + 'received at the DSN receiving station, given in elapsed '
            + 'seconds' + sd + 'after the moment specified by '
            + 'REFERENCE_TIME. Also referred to' + sd
            + 'as Earth receiving time or ERT."')
            ,
            ('"Frequency of the downlink sinusoidal' + sd + 'signal '
            + 'at the front-end of the DSN_STATION_NUMBER. It is '
            + 'computed using' + sd + 'frequency encoding information '
            + 'recorded by the RSR and the frequency' + sd
            + 'offset from the center of the recording bandwidth '
            + 'observed in' + sd + 'computed power spectra of the '
            + 'measured I/Q samples. See formula' + sd
            + 'given in Section 2.6 of JPLD-16765."')
            ,
            ('"The residual frequency is the' + sd + 'difference between '
            + 'the sky_frequency and the Doppler-shifted' + sd
            + 'frequency of the received sinusoid computed based on the'
            + sd + 'reconstructed spacecraft trajectory and smoothed '
            + 'using weighted' + sd + 'least-square spline fit across '
            + 'the extent of the ring system. Time' + sd + 'history of '
            + 'the estimate is used to steer the frequency of the' + sd
            + 'received sinusoid to the center of the recording '
            + 'bandwidth."')
            ,
            ('"Estimate of the power of the received' + sd + 'sinusoid '
            + 'in the absence of the rings. It is calculated using' + sd
            + 'weighted least-square spline fit to estimates of the '
            + 'downlink' + sd + 'carrier power in regions outside of '
            + 'Ring A and inside of Ring C' + sd + '(when feasible), '
            + 'and within prominent gaps in Ring A, the Cassini' + sd
            + 'Division, and Ring C. Time history of the estimate is '
            + 'used to' + sd + 'normalize the free-space power of the '
            + 'downlink signal to unity' + sd + 'across the extent '
            + 'of the ring system."')
            ]

    object_values = []
    object_values.append(object_names)
    object_values.append(col_num)
    object_values.append(data_types)
    object_values.append(start_bytes_list)
    object_values.append(bytes_list)
    object_values.append(formats)
    object_values.append(units)
    object_values.append(reference_times)
    object_values.append(object_descriptions)

        
    str_lbl = {
        'string_delimiter': sd,
        'alignment_column': alignment_column,
        'series_alignment_column': series_alignment_column,
        'keywords_values': keywords_values,
        'keywords_NAIF_TOOLKIT_VERSION': NAIF_TOOLKIT_VERSION_dict,
        'description': FILE_DESCRIPTION,
        'keywords_series': SERIES_dict,
        'object_keys': object_keys,
        'object_values': object_values,
        'history': HISTORY_dict
        }
    return str_lbl

def write_cal_series(rev_info, cal_inst, title, outdir, prof_dir):
    """
    This function writes a CAL series, which includes a data and label file.

    Args:
        rev_info (dict): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        cal_inst (class): Instance of Calibration class
        title (str): Name of the output .TAB and .LBL file, not including
                           extensions. Date in YYYYMMDD format will be added
                           onto series_name
        outdir (str): Path to output directory
        prof_dir (str): Direction of ring occultation for this cal_inst

    Notes:
        [1] Data entry format of %32.16F is hardcoded.
        [2] A data and label file will be output into the input "outdir"
            directory, with filenames, *YYYYMMDD.TAB and *YYYYMMDD.LBL,
            respectively, where * is "title".
    """
    current_time = time.strftime("_%Y%m%d") 
    outfile_tab = outdir + title.upper() + current_time + '.TAB'
    outfile_lbl = outdir + title.upper() + current_time + '.LBL'

    series_name = '"' + outfile_tab.split('/')[-1] + '"' 

    fmt = '%32.16F' 

    
    # Write data file
    write_cal_series_data(cal_inst, fmt, outfile_tab)

    # Get label file information
    str_lbl = get_cal_series_info(rev_info, cal_inst, series_name, prof_dir)

    # Write label file
    pds3.pds3_write_series_lbl(str_lbl, outfile_lbl)

    return None


