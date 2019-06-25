
"""
pds3_cal_series.py

Purpose: Write CAL data and label files in PDS3 format.

Dependencies:
    #. numpy
    #. time
    #. rss_ringoccs.tools.pds3_write_series_v2

Notes:
    [1] Contents of output CAL data and label files are meant to mimic
        CAL files from CORSS_8001 v2.
    
"""

from . import pds3_write_series_v2 as pds3
import numpy as np
import pdb
import time

def write_cal_series_data(cal_inst, out_file):
    """
    This writes a CAL data file with columns: observed event time,
    sky frequency, offset frequency, and fit to free-space power.

    Arguments
        :cal_inst (*class*): Instance of Calibration class
        :out_file (*str*): Path to output file
    """
    format_str = ('%14.6F,' + '%20.6F,' + '%10.6F,' + '%14.6F' + '%s')
    npts = len(cal_inst.t_oet_spm_vals)
    print('\tWriting CAL data to: \n\t\t' + out_file)
    f = open(out_file, 'w')
    for n in range(npts):
        f.write(format_str % (
            cal_inst.t_oet_spm_vals[n],
            cal_inst.f_sky_hz_vals[n],
            cal_inst.f_offset_fit_vals[n],
            cal_inst.p_free_vals[n],
            '\r\n'))

    f.close()

    return None


def get_cal_series_info(rev_info, cal_inst, series_name, prof_dir):
    """
    This returns the information needed to write a CAL label file.

    Arguments
        :rev_info (*dict*): Dictionary with keys: rsr_file, band, year, doy, 
                        dsn, occ_dir, planetary_occ_flag, rev_num
        :cal_inst (*class*): Instance of Calibration class
        :series_name (*str*): Name of the output .TAB and .LBL file, 
                            not including extensions. '_YYYYMMDD_XXXX' will
                            be added to the end of series_name 
        :prof_dir (*str*): Direction of ring occultation for this cal_inst

    Returns
        :str_lbl (*dict*): Dictionary with keys: string_delimiter,
                        alignment_column, series_alignment_column,
                        keywords_value, keywords_NAIF_TOOLKIT_VERSION,
                        description, keywords_series, object_keys,
                        object_values, history
    """
    # Get current time in _YYYYMMDD-hhmmss format
    current_time_ISOD = time.strftime("%Y-%j") + 'T' + time.strftime("%H:%M:%S")

    # Values for number of columns, number of bytes per column,
    #   number of bytes per column delimiter, number of bytes allocated to
    #   special characters per column
    formats = ['"F14.6"', '"F20.6"', '"F10.6"', '"F14.6"']
    ncol = len(formats)

    # Values for aligning equal signs
    alignment_column        = 32
    series_alignment_column = 28

    record_bytes, bytes_list, start_bytes_list = pds3.get_record_bytes(formats)

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
    geo_kernels = cal_inst.history['Positional Args']['geo_inst'][
            'Positional Args']['kernels']
    # Remove directory path from kernel list
    if isinstance(geo_kernels, list):
        geo_kernels = ['"'+x.split('/')[-1]+'"' for x in geo_kernels]
    else:
        geo_kernels = '"' + geo_kernels + '"'


    PDS_VERSION_ID = 'PDS3'
    RECORD_TYPE = 'FIXED_LENGTH'
    RECORD_BYTES = record_bytes
    FILE_RECORDS = str(len(sampling_parameter_arr))
    SERIES_NAME = series_name

    DATA_SET_ID = '"CO-SR-RSS-?/?-OCC-V0.1"'
    RING_OBSERVATION_ID = pds3.get_ring_obs_id(year, doy, band, dsn)
    PRODUCT_ID = series_name
    PRODUCT_TYPE = 'CALIBRATION_PARAMETERS'
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
    NAIF_TOOLKIT_VERSION = cal_inst.naif_toolkit_version

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
            + 'signal attributes needed to calibrate the raw ring data before'
            + sd + 'reliable optical depth profiles are computed. The '
            + 'attributes are the' + sd + 'signal sky-frequency, the '
            + 'frequency offset, and the free-space' + sd
            + 'signal power. The attributes are listed versus ''OBSERVED '
            + 'EVENT TIME''' + sd + '(Earth received time) over equal '
            + 'time increments (1 s).' + sd
            + ' ' + sd
            + 'The frequency offset estimates are required to steer the '
            + 'frequency' + sd + 'of the downlink sinusoid to a constant '
            + 'value (near the center of the' + sd 
            + 'recording bandwidth) before the sinusoid amplitude and '
            + 'phase can be' + sd + 'estimated. The free-space signal '
            + 'power estimates are required to' + sd + 'normalize the '
            + 'power of the steered sinusoid to value of about unity' + sd
            + '(optical depth of about zero) outside Ring A, inside Ring C, '
            + 'and within' + sd + 'large ring gaps. The sky-frequency '
            + 'estimates are included for completeness' + sd
            + 'and to facilitate independent checks of frequency '
            + 'calculations.' + sd
            + ' ' + sd
            + 'The radio science receiver at the DSN ground receiving station '
            + 'the (RSR)' + sd + 'steers the frequency of the received '
            + 'sinusoid so that the measured' +  sd + 'spectral line falls '
            + 'at the center of the recording bandwidth. Spectral' + sd
            + 'estimates of the measured I/Q samples can be used to '
            + 'calculate any' + sd + 'offsets of the spectral line from '
            + 'the center of the bandwidth. The' + sd
            + 'measured offsets together with other frequency steering '
            + 'information' + sd + 'encoded in the RSR recording are used '
            + 'to calculate the listed sky-' + sd
            + 'frequency based on procedures documented in the JPLD-16765. '
            + 'Part' + sd + 'of the measured frequency offset is caused '
            + 'by the use of a predicted' + sd
            + 'spacecraft trajectory to estimate '
            + 'the received Doppler-shifted signal' + sd
            + 'frequency needed to steer the RSR. This part can be removed '
            + 'by using' + sd +'the more accurate trajectory reconstructed '
            + 'by the Cassini Navigation' + sd + 'Team instead. '
            + 'The listed frequency offset is the polynomial fit '
            + 'over the' + sd + 'global extent of the rings. The free-space '
            + 'signal power estimates are' + sd + 'also calculated using '
            + 'global polynomial fits to power estimates outside' + sd
            + 'Ring A and, when possible, inside Ring C, as well as within '
            + 'major' + sd + 'ring gaps.' + sd
            + ' ' + sd
            + 'Relevant geometry calculations are based on the use of the '
            + 'Cassini' + sd + 'Navigation Team Naif Toolkit kernel files '
            + 'available at the time of' + sd + 'archiving and are listed '
            + 'above. We note that the adopted Planetary' + sd 
            + 'Constants Kernel (PCK) file is not necessarily the one '
            + 'Cassini NAV' + sd + 'associates with the listed '
            + 'reconstructed trajectory file. The' + sd 
            + 'difference this causes to estimate ring radius is well '
            + 'below 1 km' + sd + 'and has negligble impact on the '
            + 'archived products.' + sd
            + ' ' + sd
            + 'All calculations assume fixed UltraStable Oscillator (USO) '
            + 'reference' + sd + 'frequency of 8,427,222,034.34050 Hz '
            + 'at X-band, its value near the' + sd + 'beginning of the '
            + 'Cassini orbital tour. The frequency is coherently' + sd
            + 'scaled by 3/11 for S-band and by 209/55 for Ka-band. The exact '
            + 'USO' + sd + 'frequency changed slightly (at the Hz level) '
            + 'during the USO lifetime.' + sd + 'The change negligibly '
            + 'impacts the archived products. The same holds' + sd
            + 'true for the Allan deviation characterizing the stability of '
            + 'the USO.' + sd + 'Typical values of the Allan deviation '
            + 'is 2E-13 over 1 s and 1E-13' + sd + 'over 10-100 s. '
            + 'The value changed little over the USO lifetime.' + sd
            + ' ' + sd
            + 'This file was produced using the rss_ringoccs open-source '
            + 'processing suite' + sd + 'developed at Wellesley College with '
            + 'the support of the Cassini project and' + sd + 'hosted '
            + 'on GithHub at https://github.com/NASA-Planetary-Science/'
            + 'rss_ringoccs.' + sd
            + ' ' + sd
            + 'Please address any inquiries to:' + sd
            + 'Richard G. French,' + sd
            + 'Astronomy Department, Wellesley College;' + sd
            + 'Wellesley, MA 02481-8203;' + sd
            + '(781) 283-3747;' + sd
            + 'rfrench@wellesley.edu."')

    HIST_USER_NAME = cal_inst.history['User Name']
    HIST_HOST_NAME = cal_inst.history['Host Name']
    HIST_RUN_DATE = cal_inst.history['Run Date']
    HIST_PYTHON_VERSION = cal_inst.history['Python Version']
    HIST_OPERATING_SYSTEM = cal_inst.history['Operating System']
    HIST_SOURCE_DIR = cal_inst.history['Source Directory']
    HIST_SOURCE_FILE = cal_inst.history['Source File']
    HIST_INPUT_VARIABLES = cal_inst.history['Positional Args']
    HIST_INPUT_KEYWORDS = cal_inst.history['Keyword Args']
    HIST_ADD_INFO = cal_inst.history['Additional Info']
    HIST_RSSOCC_VERSION = cal_inst.history['rss_ringoccs Version']
    HIST_description = ('This is a record of the processing steps'
                        + sd + 'and inputs used to generate this file.')

    HISTORY_dict = {
            'key_order0': ['User Name', 'Host Name', 'Operating System',
                        'Python Version', 'rss_ringoccs Version']
            ,'key_order1': ['Source Directory','Source File',
                        'Positional Args', 'Keyword Args', 'Additional Info']
            , 'hist name': 'Calibration history'
            , 'User Name': HIST_USER_NAME
            , 'Host Name': HIST_HOST_NAME
            , 'Run Date': HIST_RUN_DATE
            , 'Python Version': HIST_PYTHON_VERSION
            , 'rss_ringoccs Version': HIST_RSSOCC_VERSION
            , 'Operating System': HIST_OPERATING_SYSTEM
            , 'Source Directory': HIST_SOURCE_DIR
            , 'Source File': HIST_SOURCE_FILE
            , 'Positional Args': HIST_INPUT_VARIABLES
            , 'Keyword Args': HIST_INPUT_KEYWORDS
            , 'Additional Info': HIST_ADD_INFO
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

    Arguments
        :rev_info (*dict*): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        :cal_inst (*class*): Instance of Calibration class
        :title (*str*): Name of the output .TAB and .LBL file, not including
                           extensions. Date in YYYYMMDD format and sequence
                           number in XXXX format will be added at the end
                           of series_name
        :outdir (*str*): Path to output directory
        :prof_dir (*str*): Direction of ring occultation for this cal_inst

    """
    outfile_tab = outdir + title.upper() + '.TAB'
    outfile_lbl = outdir + title.upper() + '.LBL'

    series_name = '"' + outfile_tab.split('/')[-1] + '"' 

    
    # Write data file
    write_cal_series_data(cal_inst, outfile_tab)

    # Get label file information
    str_lbl = get_cal_series_info(rev_info, cal_inst, series_name, prof_dir)

    # Write label file
    pds3.pds3_write_series_lbl(str_lbl, outfile_lbl)

    return None


