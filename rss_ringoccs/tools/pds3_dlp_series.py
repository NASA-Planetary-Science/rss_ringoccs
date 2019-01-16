
"""
pds3_dlp_series.py

Purpose: Write DLP data and label files in PDS3 format.

Dependencies:
    #. numpy
    #. time
    #. rss_ringoccs.tools.pds3_write_series_v2

Notes:
    [1] Contents of output DLP data and label files are meant to mimic DLP
        files from CORSS_8001 v2.
"""
import numpy as np
import pdb
import time
from . import pds3_write_series_v2 as pds3

def write_dlp_series_data(dlp_inst, out_file, verbose=False):
    """
    This writes a CAL data file with columns: ring radius, radius correction
    due to improved pole, radius correction due to timing offset, ring 
    longitude, observed ring azimuth, power, optical depth, phase,
    threshold optical depth, observed event time, ring event time,
    spacecraft event time, ring opening angle.

    Arguments:
        :dlp_inst (*class*): Instance of DiffractionLimitedProfile class
        :out_file (*str*): Path to output file
    """
    format_str = ('%14.6F,' + '%10.6F,' + '%10.6F,' + '%12.6F,' + '%12.6F,'
            + '%14.6E,' + '%14.6E,' + '%12.6F,' + '%14.6E,' + '%14.6F,'
            + '%14.6F,' + '%14.6F,' + '%12.6F' + '%s')


    npts = len(dlp_inst.t_oet_spm_vals)

    # Compute normalized optical depth -- NOTE: this should be added to dlp_inst
    #   as an attribute
    tau_norm_vals = -np.sin(abs(dlp_inst.B_rad_vals)) * np.log(
            dlp_inst.p_norm_vals)
    if verbose:
        print('\nWriting DLP data to: ', out_file, '\n')

    f = open(out_file, 'w')
    for n in range(npts):
        f.write(format_str % (
            dlp_inst.rho_km_vals[n],
            dlp_inst.rho_corr_pole_km_vals[n],
            dlp_inst.rho_corr_timing_km_vals[n],
            np.degrees(dlp_inst.phi_rl_rad_vals[n]),
            np.degrees(dlp_inst.phi_rad_vals[n]),
            dlp_inst.p_norm_vals[n],
            tau_norm_vals[n],
            np.degrees(dlp_inst.phase_rad_vals[n]),
            dlp_inst.tau_threshold_vals[n],
            dlp_inst.t_oet_spm_vals[n],
            dlp_inst.t_ret_spm_vals[n],
            dlp_inst.t_set_spm_vals[n],
            np.degrees(dlp_inst.B_rad_vals[n]),
            '\r\n'))
            
    f.close()


    return None

def get_dlp_series_info(rev_info, dlp_inst, series_name, prof_dir):
    """
    This returns the information needed to write a DLP label file.

    Arguments:
        :rev_info (*dict*): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        :dlp_inst (*class*): Instance of DiffractionLimitedProfile
        :series_name (*str*): Name of the output .TAB and .LBL file, 
                            not including extensions. Date in YYYYMMDD format
                            and sequence number in XXXX format will be 
                            added onto series_name
        :prof_dir (*str*): Direction of ring occultation for this dlp_inst

    Returns:
        :str_lbl (*dict*): Dictionary with keys: string_delimiter,
                        alignment_column, series_alignment_column,
                        keywords_value, keywords_NAIF_TOOLKIT_VERSION,
                        description, keywords_series, object_keys,
                        object_values, history
    """
    current_time_ISOD = time.strftime("%Y-%j") + 'T' + time.strftime("%H:%M:%S")

    # Values for number of columns, number of bytes per column,
    #   number of bytes per column delimiter, number of bytes allocated to
    #   special characters per column
    formats = ['"F14.6"', '"F10.6"', '"F10.6"', '"F12.6"', '"F12.6"'
            , '"E14.6"' , '"E14.6"' , '"F12.6"' , '"E14.6"' , '"F14.6"'
            , '"F14.6"' , '"F14.6"' , '"F12.6"']
    ncol = len(formats)

    # Values for aligning equal signs
    alignment_column        = 37
    series_alignment_column = 28

    record_bytes, bytes_list, start_bytes_list = pds3.get_record_bytes(formats)

    col_num = list(range(1, ncol+1))

    rsr_file = rev_info['rsr_file']
    band = rev_info['band']
    year = rev_info['year']
    doy = rev_info['doy']
    dsn = rev_info['dsn']

    # Extract relevant information from dlp instance
    sampling_parameter_arr = dlp_inst.rho_km_vals
    t_oet_spm_start = dlp_inst.t_oet_spm_vals[0]
    t_oet_spm_end = dlp_inst.t_oet_spm_vals[-1]

    t_ret_spm_start = dlp_inst.t_ret_spm_vals[0]
    t_ret_spm_end = dlp_inst.t_ret_spm_vals[-1]

    t_set_spm_start = dlp_inst.t_set_spm_vals[0]
    t_set_spm_end = dlp_inst.t_set_spm_vals[-1]


    PDS_VERSION_ID = 'PDS3'
    RECORD_TYPE = 'FIXED_LENGTH'
    RECORD_BYTES = record_bytes
    FILE_RECORDS = str(len(sampling_parameter_arr))
    SERIES_NAME = series_name

    DATA_SET_ID = '"CO-SR-RSS-?/?-OCC-V0.1"'
    RING_OBSERVATION_ID = pds3.get_ring_obs_id(year, doy, band, dsn)
    PRODUCT_ID = series_name
    PRODUCT_TYPE = 'RING_PROFILE'
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
    RING_EVENT_START_TIME = pds3.get_ISOD_str(t_ret_spm_start, year, doy)
    RING_EVENT_STOP_TIME = pds3.get_ISOD_str(t_ret_spm_end, year, doy)
    SPACECRAFT_EVENT_START_TIME = pds3.get_ISOD_str(t_set_spm_start, year, doy)
    SPACECRAFT_EVENT_STOP_TIME = pds3.get_ISOD_str(t_set_spm_end, year, doy)
    SPACECRAFT_CLOCK_START_COUNT = '"UNK"'
    SPACECRAFT_CLOCK_STOP_COUNT = '"UNK"'
    REVOLUTION_NUMBER = rev_info['rev_num']
    DSN_STATION_NUMBER = dsn.split('-')[-1]

    OBSERVATION_TYPE = '"OCCULTATION PROFILE"'
    OCCULTATION_TYPE = 'RADIO'
    FEATURE_NAME = '"RING SYSTEM"'
    PLANETARY_OCCULTATION_FLAG = rev_info['planetary_occ_flag']
    RING_OCCULTATION_DIRECTION = rev_info['occ_dir']
    RING_PROFILE_DIRECTION = prof_dir
    FREQUENCY_BAND = band

    wavelength_dict = {
              '"X"': "3.5574   <cm>"
            , '"K"': "0.93617   <cm>"
            , '"S"': "13.044   <cm>"
            }

    WAVELENGTH = wavelength_dict[band]
    RADIAL_RESOLUTION = str(float(pds3.get_sampling_interval(
        sampling_parameter_arr))*2.) + '   <km>'
    RADIAL_SAMPLING_INTERVAL = pds3.get_sampling_interval(
            sampling_parameter_arr) + '   <km>'
    MINIMUM_RING_RADIUS = str(
            round(min(dlp_inst.rho_km_vals), 4)) + '   <km>'
    MAXIMUM_RING_RADIUS = str(
            round(max(dlp_inst.rho_km_vals),4)) + '   <km>'
    MINIMUM_RING_LONGITUDE = str(
            round(min(np.degrees(dlp_inst.phi_rl_rad_vals)), 4)) + '   <deg>'
    MAXIMUM_RING_LONGITUDE = str(
            round(max(np.degrees(dlp_inst.phi_rl_rad_vals)), 4)) + '   <deg>'
    MINIMUM_OBSERVED_RING_AZIMUTH = str(
            round(min(np.degrees(dlp_inst.phi_rad_vals)), 4)) + '   <deg>'
    MAXIMUM_OBSERVED_RING_AZIMUTH = str(
            round(max(np.degrees(dlp_inst.phi_rad_vals)), 4)) + '   <deg>'
    MINIMUM_OBSERVED_RING_ELEVATION = str(
            round(min(np.degrees(dlp_inst.B_rad_vals)), 4)) + '   <deg>'
    MAXIMUM_OBSERVED_RING_ELEVATION = str(
            round(max(np.degrees(dlp_inst.B_rad_vals)), 4)) + '   <deg>'
    LOWEST_DETECTABLE_OPACITY = str(min(dlp_inst.tau_threshold_vals))
    HIGHEST_DETECTABLE_OPACITY = str(max(dlp_inst.tau_threshold_vals))


    NAIF_TOOLKIT_VERSION = ''
    SPICE_FILE_NAME = ''


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
                         , 'RING_EVENT_START_TIME'
                         , 'RING_EVENT_STOP_TIME'
                         , 'SPACECRAFT_EVENT_START_TIME'
                         , 'SPACECRAFT_EVENT_STOP_TIME'
                         , 'SPACECRAFT_CLOCK_START_COUNT'
                         , 'SPACECRAFT_CLOCK_STOP_COUNT'
                         , 'REVOLUTION_NUMBER'
                         , 'DSN_STATION_NUMBER']
            , 'INSTRUMENT_HOST_NAME': INSTRUMENT_HOST_NAME
            , 'INSTRUMENT_HOST_ID': INSTRUMENT_HOST_ID
            , 'INSTRUMENT_NAME' : INSTRUMENT_NAME
            , 'INSTRUMENT_ID': INSTRUMENT_ID
            , 'MISSION_PHASE_NAME': MISSION_PHASE_NAME
            , 'TARGET_NAME': TARGET_NAME
            , 'START_TIME' : START_TIME
            , 'STOP_TIME' : STOP_TIME
            , 'RING_EVENT_START_TIME': RING_EVENT_START_TIME
            , 'RING_EVENT_STOP_TIME': RING_EVENT_STOP_TIME
            , 'SPACECRAFT_EVENT_START_TIME': SPACECRAFT_EVENT_START_TIME
            , 'SPACECRAFT_EVENT_STOP_TIME': SPACECRAFT_EVENT_STOP_TIME
            , 'SPACECRAFT_CLOCK_START_COUNT': SPACECRAFT_CLOCK_START_COUNT
            , 'SPACECRAFT_CLOCK_STOP_COUNT': SPACECRAFT_CLOCK_STOP_COUNT
            , 'REVOLUTION_NUMBER': REVOLUTION_NUMBER
            , 'DSN_STATION_NUMBER': DSN_STATION_NUMBER
            }

    NAIF_TOOLKIT_VERSION_dict = {
            'key_order': ['NAIF_TOOLKIT_VERSION','SPICE_FILE_NAME']
            , 'NAIF_TOOLKIT_VERSION': NAIF_TOOLKIT_VERSION
            , 'SPICE_FILE_NAME': SPICE_FILE_NAME
            }

    OBSERVATION_TYPE_dict = {
             'key_order': ['OBSERVATION_TYPE', 'OCCULTATION_TYPE',
                        'FEATURE_NAME', 'PLANETARY_OCCULTATION_FLAG',
                        'RING_OCCULTATION_DIRECTION', 'RING_PROFILE_DIRECTION',
                        'FREQUENCY_BAND', 'WAVELENGTH']
             , 'OBSERVATION_TYPE': OBSERVATION_TYPE
             , 'OCCULTATION_TYPE': OCCULTATION_TYPE
             , 'FEATURE_NAME': FEATURE_NAME
             , 'PLANETARY_OCCULTATION_FLAG': PLANETARY_OCCULTATION_FLAG
             , 'RING_OCCULTATION_DIRECTION': RING_OCCULTATION_DIRECTION
             , 'RING_PROFILE_DIRECTION': RING_PROFILE_DIRECTION
             , 'FREQUENCY_BAND': FREQUENCY_BAND
             , 'WAVELENGTH': WAVELENGTH
             }

    RADIAL_RESOLUTION_dict = {
             'key_order': ['RADIAL_RESOLUTION', 'RADIAL_SAMPLING_INTERVAL',
                        'MINIMUM_RING_RADIUS', 'MAXIMUM_RING_RADIUS',
                        'MINIMUM_RING_LONGITUDE', 'MAXIMUM_RING_LONGITUDE',
                        'MINIMUM_OBSERVED_RING_AZIMUTH',
                        'MAXIMUM_OBSERVED_RING_AZIMUTH',
                        'MINIMUM_OBSERVED_RING_ELEVATION',
                        'MAXIMUM_OBSERVED_RING_ELEVATION',
                        'LOWEST_DETECTABLE_OPACITY',
                        'HIGHEST_DETECTABLE_OPACITY']
             , 'RADIAL_RESOLUTION': RADIAL_RESOLUTION
             , 'RADIAL_SAMPLING_INTERVAL': RADIAL_SAMPLING_INTERVAL
             , 'MINIMUM_RING_RADIUS': MINIMUM_RING_RADIUS
             , 'MAXIMUM_RING_RADIUS': MAXIMUM_RING_RADIUS
             , 'MINIMUM_RING_LONGITUDE': MINIMUM_RING_LONGITUDE
             , 'MAXIMUM_RING_LONGITUDE': MAXIMUM_RING_LONGITUDE
             , 'MINIMUM_OBSERVED_RING_AZIMUTH': MINIMUM_OBSERVED_RING_AZIMUTH
             , 'MAXIMUM_OBSERVED_RING_AZIMUTH': MAXIMUM_OBSERVED_RING_AZIMUTH
             , 'MINIMUM_OBSERVED_RING_ELEVATION': 
                        MINIMUM_OBSERVED_RING_ELEVATION
             , 'MAXIMUM_OBSERVED_RING_ELEVATION': 
                        MAXIMUM_OBSERVED_RING_ELEVATION
             , 'LOWEST_DETECTABLE_OPACITY': LOWEST_DETECTABLE_OPACITY
             , 'HIGHEST_DETECTABLE_OPACITY': HIGHEST_DETECTABLE_OPACITY
             }


    qq = "'"

    sd = '|'
    FILE_DESCRIPTION = ('"The ''DLP'' file contains' + sd
            + 'calibrated but diffraction-limited optical depth and phase '
            + 'shift profiles' + sd + 'of Saturn''s rings, that is, '
            + 'calibrated profiles before reconstruction' + sd
            + 'to remove diffraction effects. The frequency/phase '
            + 'measurements' + sd + 'reference is the Cassini UltraStable '
            + 'Oscillator (USO). The diffraction' + sd + 'reconstruction '
            + 'is carried out using algorithms described at length in' + sd
            + 'MAROUFETAL1986. Practical implementation steps are provided '
            + 'in an' + sd + 'included documentation file.' + sd
            + ' ' + sd
            + 'Several additional companion products share the same '
            + 'RING_OBSERVATION_ID' + sd + '(listed in the header of the LBL '
            + 'file) as this product. These include' + sd + 'one or more '
            + 'reconstructed (TAU) files for the same occultation but' + sd
            + 'at different resolutions, and two which provide geometry and '
            + 'calibration' + sd + 'data. The latter two have file names '
            + 'constructed from the same root as' + sd
            + 'this file with the field for radial resolution removed '
            + 'and the ''DLP''' + sd + 'replaced by either ''GEO'' '
            + '(geometry) or ''CAL'' (carrier frequency and' + sd
            + 'power calibration).' + sd
            + ' ' + sd
            + 'We define spatial resolution as the shortest resolvable '
            + 'wavelength in the' + sd + 'DLP profiles. It''s the inverse '
            + 'of the highest spatial frequency preserved' + sd
            + 'in the data, which is determined by the bandwidth of '
            + 'the lowpass filter' + sd + 'used to decimate the calibrated '
            + 'data. The archived DLP files have a shortest' + sd
            + 'resolvable wavelength of 0.5 km (the ''500M'' descriptor '
            + 'in the file name),' + sd + 'and are sampled every 0.25 km '
            + '(sampled at the Nyquist rate). They are' + sd
            + 'provided for the benefit of users wanting to experiment with '
            + 'their own' + sd + 'diffraction reconstruction implementation. '
            + 'Recovered profiles should look' + sd + 'identical to the '
            + 'archived 1 km resolution reconstructed profiles provided' + sd
            + 'that the reconstruction is implemented to achieve'
            + ' ''processing resolution'' of' + sd + '0.75 km, where '
            + 'the latter is defined by Eq. 19 of MAROUFETAL1986.' + sd
            + ' ' + sd
            + 'All archived data products were generated assuming '
            + 'fixed USO '
            + 'reference' + sd + 'frequency of 8,427,222,034.34050 Hz '
            + 'at X-band, its value near the beginning' + sd + 'of the '
            + 'Cassini orbital tour. The frequency is coherently '
            + 'scaled by 3/11' + sd + 'for S-band and by 209/55 for Ka-band. '
            + 'The exact USO frequency changed' + sd
            + 'slightly (at the Hz level) during the USO lifetime. '
            + 'The change negligibly ' + sd
            + 'impacts the archived products. The same holds '
            + 'true for the Allan deviation' + sd + 'characterizing the '
            + 'stability of the USO. Typical values of the Cassini USO' + sd
            + 'Allan deviation is 2E-13 over 1 s and 1E-13 over 10-100 s. '
            + 'The USO Allan' + sd + 'deviation changed little over the '
            + 'USO lifetime.' + sd
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
            
    HIST_USER_NAME = dlp_inst.history['User Name']
    HIST_HOST_NAME = dlp_inst.history['Host Name']
    HIST_RUN_DATE = dlp_inst.history['Run Date']
    HIST_PYTHON_VERSION = dlp_inst.history['Python Version']
    HIST_OPERATING_SYSTEM = dlp_inst.history['Operating System']
    HIST_SOURCE_DIR = dlp_inst.history['Source Directory']
    HIST_SOURCE_FILE = dlp_inst.history['Source File']
    HIST_INPUT_VARIABLES = dlp_inst.history['Positional Args']
    HIST_INPUT_KEYWORDS = dlp_inst.history['Keyword Args']
    HIST_ADD_INFO = dlp_inst.history['Additional Info']
    HIST_RSSOCC_VERSION = dlp_inst.history['rss_ringoccs Version']
    HIST_description = ('This is a record of the processing steps'
                        + sd + 'and inputs used to generate this file.')

    HISTORY_dict = {
            'key_order0': ['User Name', 'Host Name', 'Operating System',
                        'Python Version', 'rss_ringoccs Version']
            ,'key_order1': ['Source Directory','Source File',
                        'Positional Args', 'Keyword Args', 'Additional Info']
            , 'hist name': 'DiffractionLimitedProfile history'
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
        blank,
        OBSERVATION_TYPE_dict,
        blank,
        RADIAL_RESOLUTION_dict,
        blank
        ]

    SERIES = 'SERIES'
    SERIES_NAME = '"OCCULTATION PROFILE"'
    SERIES_INTERCHANGE_FORMAT = 'ASCII'
    SERIES_COLUMNS = '12'
    SERIES_ROWS = FILE_RECORDS
    SERIES_ROW_BYTES = RECORD_BYTES
    SERIES_SAMPLING_PARAMETER_NAME = '"RING RADIUS"'
    SERIES_SAMPLING_PARAMETER_UNIT = '"KILOMETER"'
    SERIES_MINIMUM_SAMPLING_PARAMETER = str(min(sampling_parameter_arr))
    SERIES_MAXIMUM_SAMPLING_PARAMETER = str(max(sampling_parameter_arr))
    SERIES_SAMPLING_PARAMETER_INTERVAL = pds3.get_sampling_interval(
            sampling_parameter_arr)
    SERIES_DESCRIPTION = ('"This series contains variations of the' + sd
            + 'optical depth and phase shift profiles as a function of '
            + 'RING RADIUS. Also' + sd + 'included are'
            + ' some pertinent event times, and geometry parameters."')

    
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
            '"RING RADIUS"'
            , '"RADIUS CORRECTION DUE TO IMPROVED POLE"'
            , '"RADIUS CORRECTION DUE TO TIMING OFFSET"'
            , '"RING LONGITUDE"'
            , '"OBSERVED_RING_AZIMUTH"'
            , '"NORMALIZED SIGNAL POWER"'
            , '"NORMAL OPTICAL DEPTH"'
            , '"PHASE SHIFT"'
            , '"NORMAL OPTICAL DEPTH THRESHOLD"'
            , '"OBSERVED EVENT TIME"'
            , '"RING EVENT TIME"'
            , '"SPACECRAFT EVENT TIME"'
            , '"OBSERVED RING ELEVATION"'
            ]
    n_objects = len(object_names)
    data_types = ['ASCII_REAL'] * n_objects
    formats = ['"F32.16"'] * n_objects
    units = ['"KILOMETER"', '"N/A"', '"N/A"', '"DEGREE"',
            '"DEGREE"', '"N/A"', '"N/A"', '"DEGREE"', '"N/A"', '"SECOND"',
            '"SECOND"', '"SECOND"', '"DEGREE"']

    es = ''

    reference_times = [es, es, es, es, es, es, es, es, es,
            OBJECT_REFERENCE_TIME, OBJECT_REFERENCE_TIME,
            OBJECT_REFERENCE_TIME,es]

    object_descriptions = [
            ('"Radial distance from the center of' + sd + 'Saturn '
            + 'to the ring-plane intercept point at the RING EVENT TIME."')
            ,
            ('"This is a placeholder column to match' + sd
            + 'the ''DLP'' file format provided in: Cassini_RSS_Ring_'
            + 'Profiles_2018_Archive.' + sd + 'There is no information '
            + 'in this column -- it is populated with zeros."')
            ,
            ('"This is a placeholder column to match' + sd
            + 'the ''DLP'' file format provided in: Cassini_RSS_Ring_'
            + 'Profiles_2018_Archive.' + sd + 'There is no information '
            + 'in this column -- it is populated with zeros."')
            ,
            ('"Inertial (J2000) longitude in the ring' + sd 
            + 'plane of the ring-plane intercept point at the '
            + 'RING EVENT TIME."')
            ,
            ('"Measured at the ring-plane intercept' + sd + 'point, '
            + 'starting from the direction of a photo heading to the' + sd
            + 'observer (Earth receiving station), and ending at the '
            + 'direction of' + sd + 'a local radial vector. This angle '
            + 'is projected into the ring plane' + sd + 'and measured '
            + 'in the prograde direction. As seen from the observer,' + sd
            + 'it equals 90 degrees along the ring ansa and 270 degrees '
            + 'along the' + sd + 'left ansa. Values range from 0 to 360 '
            + 'degrees. This convention for' + sd + 'observed ring '
            + 'azimuth differs from that adopted in MAROUFETAL1986 by' + sd
            + '180 degrees."')
            ,
            ('"Power (amplitude square) of the measured' + sd
                + 'ring-attenuated diffraction-limited radio signal, '
                + 'normalized by its' + sd + 'value in the absence of the '
                + 'rings (normalized to unity in free-space).' + sd
                + 'The value may be used to compute the measured '
                + 'diffraction-limited' + sd + 'oblique optical depth '
                + 'as the negative natural logarithm of the NORMALIZED' + sd
                + 'SIGNAL POWER."')
            ,
            ('"The normal optical depth obtained' + sd + 'from its measured '
                + 'oblique value csaled by the sine of the absolute' + sd
                + 'value of ring opening angle (OBSERVED RING ELEVATION)')
            ,
            ('"The difference between the phase of' + sd + 'the coherent '
            + 'sinusoid passing directly through the rings (the direct'
            + sd + 'signal) and its value had the rings been absent. The '
            + 'phase shift is' + sd + 'computed from the measured '
            + 'complex ring transmittance;' + sd
            + 'see Eq. (22) of MAROUFETAL1986. The phase reference '
            + 'of the original' + sd + 'measurements is the Cassini '
            + 'UltraStable Oscillator, or USO."')
            ,
            ('"The value of the normal optical depth' + sd + 'for which '
            + 'the magnitude of the measured complex ring transmittance'
            + sd + 'is equal in numerical value to the 70% confidence '
            + 'interval of the' + sd + 'noisy measurement (signal-to-'
            + 'noise ratio, or SNR, of about unity).' + sd
            + 'See Eq. (26) of MAROUFETAL1986. The threshold value '
            + 'and corresponding' + sd + 'LOWEST_DETECTABLE_OPACITY '
            + 'and HIGHEST_DETECTABLE_OPACITY are based on' + sd
            + 'the measurement SNR as defined in MAROUFETAL, Eq. (27), '
            + 'and include' + sd + 'only the contribution of thermal '
            + 'noise. Further limitations on the' + sd + 'lowest '
            + 'detectable optical depth may be imposed by ground antenna'
            + sd + 'pointing errors. The imposed limitations may be '
            + 'assessed from' + sd + ' examination of the optical depth '
            + 'fluctuations observed during the' + sd + 'free-space '
            + 'baseline period over the radial scales of interest."')
            ,
            ('"The instant at which photons were' + sd + 'received at '
            + 'the DSN receiving station, given in elapsed seconds' + sd
            + 'after the moment specified by REFERENCE_TIME. Also '
            + 'referred to' + sd + 'as Earth receiving time or ERT."')
            ,
            ('"The time at which photos left the' + sd + 'ring plane. '
            + 'This time is earlier than the associated' + sd
            + 'OBSERVED EVENT TIME by an amount equal to the light travel '
            + 'time.' + sd + 'RING EVENT TIME is given in elapsed seconds '
            + 'after the moment' + sd + 'specified by REFERENCE_TIME."')
            ,
            ('"The time at which photons left the' + sd + 'spacecraft, '
            + 'given in elapsed seconds after the moment specified by' 
            + sd + 'REFERENCE_TIME. Also referred to as SCET."')
            ,
            ('"The angle measured at the ring' + sd + 'intercept point, '
            + 'starting from the ring plane and ending in the' + sd
            + 'direction of the photon heading toward the observer. '
            + 'This angle is' + sd + 'positive on the north side of '
            + 'Saturn''s rings and negative on the' + sd + 'south side. '
            + 'Its value is nearly constant over the duration of a ' + sd
            + 'ring occultation experiment and is nearly equal to the '
            + 'ring opening' + sd + 'angle (Earth elevation angle '
            + 'above the ring plane)."')
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


def write_dlp_series(rev_info, dlp_inst, title, outdir, prof_dir):
    """
    This function writes a DLP series, which includes a data and label file.

    Arguments:
        :rev_info (*dict*): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        :dlp_inst (*class*): Instance of DiffractionLimitedProfile class
        :title (*str*): Name of the output .TAB and .LBL file, not including
                           extensions. Date in YYYYMMDD format and sequence
                           number in XXXX format will be added to series_name
        :outdir (*str*): Path to output directory
        :prof_dir (*str*): Direction of ring occultation for this dlp_inst

    """


    outfile_tab = outdir + title.upper() + '.TAB'
    outfile_lbl = outdir + title.upper() + '.LBL'

    series_name = '"' + outfile_tab.split('/')[-1] + '"'

    # Write data file
    write_dlp_series_data(dlp_inst, outfile_tab)

    # Get label file information
    str_lbl = get_dlp_series_info(rev_info, dlp_inst, series_name, prof_dir)

    # Write label file
    pds3.pds3_write_series_lbl(str_lbl, outfile_lbl)
    return None
