
'''
pds3_geo_series.py

Purpose: Write GEO data and label files in PDS3 format.

Dependencies:
    #. numpy
    #. time
    #. rss_ringoccs.tools.pds3_write_series_v2

Notes:
    [1] Contents of output GEO data and label files are meant to mimic
        GEO files from CORSS_8001 v2.
'''
import pdb
import time
from . import pds3_write_series_v2 as pds3
import numpy as np


def write_geo_series_data(geo_inst, out_file):
    """
    This writes a GEO data file with columns: observed event time, ring
    event time, spacecraft event time, ring radius, ring longitude,
    observed ring azimuth, ring opening angle, distance from spacecraft to
    ring intercept point, radial velocity, azimuthal velocity, Fresnel scale,
    impact radius, x-component of spacecraft position, y-component of
    spacecraft position, z-component of spacecraft position, x-component
    of spacecraft velocity, y-component of spacecraft velocity, z-component
    of spacecraft velocity, observed spacecraft latitude.

    Arguments:
        :geo_inst (*class*): Instance of Geometry class
        :out_file (*str*): Path to output file
    """
    format_str = ('%14.6F,'*4 + '%12.6F,'*3 + '%16.6F,' + '%14.6F,'*4
                  + '%16.6F,'*3 + '%14.6F,'*3 + '%12.6F' + '%s')

    npts = len(geo_inst.t_oet_spm_vals)

    print('\tWriting GEO data to: \n\t\t' + out_file)
    f = open(out_file, 'w')

    for n in range(npts):
        f.write(format_str % (
            geo_inst.t_oet_spm_vals[n],
            geo_inst.t_ret_spm_vals[n],
            geo_inst.t_set_spm_vals[n],
            geo_inst.rho_km_vals[n],
            geo_inst.phi_rl_deg_vals[n],
            geo_inst.phi_ora_deg_vals[n],
            geo_inst.B_deg_vals[n],
            geo_inst.D_km_vals[n],
            geo_inst.rho_dot_kms_vals[n],
            geo_inst.phi_rl_dot_kms_vals[n],
            geo_inst.F_km_vals[n],
            geo_inst.R_imp_km_vals[n],
            geo_inst.rx_km_vals[n],
            geo_inst.ry_km_vals[n],
            geo_inst.rz_km_vals[n],
            geo_inst.vx_kms_vals[n],
            geo_inst.vy_kms_vals[n],
            geo_inst.vz_kms_vals[n],
            geo_inst.elev_deg_vals[n],
            '\r\n'))
    f.close()
    return None


def get_geo_series_info(rev_info, geo_inst, series_name, prof_dir):
    """
    This returns the information needed to write a GEO label file.

    Arguments
        :rev_info (*dict*): Dictionary with keys: rsr_file, band, year, doy,
                        dsn, occ_dir, planetary_occ_flag, rev_num
        :geo_inst (*class*): Instance of Geometry class
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
    # Get current time in ISOD format
    current_time_ISOD = time.strftime("%Y-%j") + 'T' + time.strftime("%H:%M:%S")

    # Values for number of columns, number of bytes per column,
    #   number of bytes per column delimiter, number of bytes allocated to
    #   special characters per column
    formats = ['"F14.6"', '"F14.6"','"F14.6"', '"F14.6"',
            '"F12.6"', '"F12.6"', '"F12.6"',
            '"F16.6"',
            '"F14.6"', '"F14.6"', '"F14.6"', '"F14.6"',
            '"F16.6"', '"F16.6"', '"F16.6"',
            '"F14.6"', '"F14.6"', '"F14.6"',
            '"F12.6"']
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
    sampling_parameter_arr = geo_inst.t_oet_spm_vals
    t_oet_spm_start = geo_inst.t_oet_spm_vals[0]
    t_oet_spm_end = geo_inst.t_oet_spm_vals[-1]
    geo_kernels = geo_inst.kernels
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
    PRODUCT_TYPE = 'OCCULTATION_GEOMETRY'
    PRODUCT_CREATION_TIME = current_time_ISOD
    PRODUCER_ID = '"TC2017"'

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

    NAIF_TOOLKIT_VERSION = geo_inst.naif_toolkit_version
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
                        'PRODUCT_CREATION_TIME', 'PRODUCER_ID']
            , 'DATA_SET_ID': DATA_SET_ID
            , 'RING_OBSERVATION_ID': RING_OBSERVATION_ID
            , 'PRODUCT_ID': PRODUCT_ID
            , 'PRODUCT_TYPE': PRODUCT_TYPE
            , 'PRODUCT_CREATION_TIME': PRODUCT_CREATION_TIME
            , 'PRODUCER_ID': PRODUCER_ID 
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

    # Description of file, with sd marking end of lines
    FILE_DESCRIPTION = ('"This file contains parameters' + sd
            + 'pertinent to the ring occultation observation '
            + 'geometry needed for' + sd + 'reconstruction of '
            + 'the diffraction limited profiles and for other' + sd
            + 'data analysis purposes. Light-time effects have been '
            + 'accounted for' + sd + 'in computing these parameters. '
            + 'For more details, see MAROUFETAL1986' + sd
            + 'and Appendix B of ROSEN1989.' + sd 
            + ' ' + sd
            + 'The geometry calculations are based on the use of the official'
            + ' Cassini' + sd + 'Navigation Team NAIF Toolkit kernel '
            + 'files available at the time of' + sd
            + 'archiving and are listed above. We note that the adopted '
            + 'Planetary' + sd + 'Constants Kernel (PCK) file is not '
            + 'necessarily the one Cassini NAV' + sd
            + 'associates with the listed reconstructed trajectory file. The'
            + sd + 'difference this causes to estimated ring radius is well '
            + 'below 1 km and' + sd + 'has negligible impact on the archived '
            + 'products.' + sd
            + ' ' + sd
            + 'All calculations assumed fixed UltraStable Oscillator (USO) '
            + 'reference' + sd + 'frequency of 8,427,222,034.34050 Hz at '
            + 'X-band, its value near the' + sd
            + 'beginning of the Cassini orbital tour. The frequency is '
            + 'coherently' + sd + 'scaled by 3/11 for S-band and 209/55 '
            + 'for Ka-band. The exact USO' + sd 
            + 'requency changed slightly (at the Hz level) during the USO '
            + 'lifetime.' + sd + 'The change negligibly impacts the '
            + 'archived products. The same holds' + sd
            + 'true for the Allan deviation characterizing the stability '
            + 'of the USO.' + sd + 'Typical values of the Allan deviation is '
            + '2E-13 over 1 s and 1E-13 over' + sd + '10-100 s. '
            + 'The values changed little over the lifetime of the USO.' + sd
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

    # Extract history information
    HIST_USER_NAME = geo_inst.history['User Name']
    HIST_HOST_NAME = geo_inst.history['Host Name']
    HIST_RUN_DATE = geo_inst.history['Run Date']
    HIST_PYTHON_VERSION = geo_inst.history['Python Version']
    HIST_OPERATING_SYSTEM = geo_inst.history['Operating System']
    HIST_SOURCE_DIR = geo_inst.history['Source Directory']
    HIST_SOURCE_FILE = geo_inst.history['Source File']
    HIST_INPUT_VARIABLES = geo_inst.history['Positional Args']
    HIST_INPUT_KEYWORDS = geo_inst.history['Keyword Args']
    HIST_ADD_INFO = geo_inst.history['Additional Info']
    HIST_RSSOCC_VERSION = geo_inst.history['rss_ringoccs Version']
    HIST_description = ('This is a record of the processing steps'
                        + sd + 'and inputs used to generate this file.')

    HISTORY_dict = {
            'key_order0': ['User Name', 'Host Name', 'Operating System',
                        'Python Version', 'rss_ringoccs Version']
            ,'key_order1': ['Source Directory','Source File',
                        'Positional Args', 'Keyword Args', 'Additional Info']
            , 'hist name': 'Geometry history'
            , 'User Name': HIST_USER_NAME
            , 'Host Name': HIST_HOST_NAME
            #, 'Run Date': HIST_RUN_DATE
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


    # Write series information
    SERIES = 'SERIES'
    SERIES_NAME = '"OCCULTATION GEOMETRY"'
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
    SERIES_DESCRIPTION = ('"Series contains ring occultation' + sd
            + 'geometry data needed for diffraction reconstruction '
            + 'of the rings' + sd + 'optical depth profile and for '
            + 'other data analysis purposes."')

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
            , '"RING EVENT TIME"'
            , '"SPACECRAFT EVENT TIME"'
            , '"RING RADIUS"'
            , '"RING LONGITUDE"'
            , '"OBSERVED RING AZIMUTH"'
            , '"OBSERVED RING ELEVATION"'
            , '"SPACECRAFT TO RING INTERCEPT DISTANCE"'
            , '"RING INTERCEPT RADIAL VELOCITY"'
            , '"RING INTERCEPT AZIMUTHAL VELOCITY"'
            , '"FRESNEL SCALE"'
            , '"IMPACT RADIUS"'
            , '"SPACECRAFT POSITION X"'
            , '"SPACECRAFT POSITION Y"'
            , '"SPACECRAFT POSITION Z"'
            , '"SPACECRAFT VELOCITY X"'
            , '"SPACECRAFT VELOCITY Y"'
            , '"SPACECRAFT VELOCITY Z"'
            , '"OBSERVED SPACECRAFT LATITUDE"'
            ]

    n_objects = len(object_names)
    data_types = ['ASCII_REAL'] * n_objects
    units = ['"SECOND"', '"SECOND"', '"SECOND"',
            '"KILOMETER"', '"DEGREE"', '"DEGREE"', '"DEGREE"',
            '"KILOMETER"', '"KM/SEC"', '"KM/SEC"',
            '"KILOMETER"', '"KILOMETER"', '"KILOMETER"', '"KILOMETER"',
            '"KILOMETER"', '"KM/SEC"', '"KM/SEC"', '"KM/SEC"',
            '"DEGREE"']
    
    es = ''
    reference_times = [OBJECT_REFERENCE_TIME, OBJECT_REFERENCE_TIME,
            OBJECT_REFERENCE_TIME, es, es, es, es, es, es, es, es, es, 
            es, es, es, es, es, es, es]
    
    object_descriptions = [
            ('"The instant at which photons were' + sd
            + 'received at the DSN receiving station, given in '
            + 'elapsed seconds' + sd + 'after the moment specified '
            + 'by REFERENCE_TIME. Also referred to' + sd
            + 'as Earth receiving time or ERT."')
            ,
            ('"The time at which photons left the' + sd + 'ring plane. '
            + 'This time is earlier than the associated' + sd
            + 'OBSERVED EVENT TIME by an amount equal to the '
            + 'light travel time.' + sd + 'RING EVENT TIME is given '
            + 'in elapsed seconds after the moment' + sd
            + 'specified by REFERENCE_TIME"')
            ,
            ('"The time at which photons left the' + sd
            + 'spacecraft, given in elapsed seconds after '
            + 'the moment specified by' + sd + 'REFERENCE_TIME. '
            + 'Also referred to as SCET."')
            ,
            ('"Radial distance from the center of' + sd
            + 'Saturn to the ring-plane intercept point at '
            + 'the RING EVENT TIME."')
            ,
            ('"Inertial (J2000) longitude in the ring' + sd + 'plane '
            + 'of the ring-plane intercept point at the RING EVENT TIME."')
            ,
            ('"Measured at the ring-plane intercept' + sd 
            + 'point, starting from the direction of a photon heading to '
            + 'the' + sd + 'observer (Earth receiving station), '
            + 'and ending at the direction of' + sd
            + 'a local radial vector. This angle is '
            + 'projected into the ring plane' + sd + 'and measured in the '
            + 'prograde direction. As seen from the observer, ' + sd
            + 'it equals 90 degrees along the right ansa and 270 degrees '
            + 'along the' + sd + 'left ansa. Values range from 0 '
            + 'to 360 degrees. This convention for' + sd
            + 'observed ring azimuth differs from that adopted in '
            + 'MAROUFETAL1986 by' + sd + '180 degrees."')
            ,
            ('"The angle measured at the ring' + sd
            + 'intercept point, starting from the ring plane '
            + 'and ending in the' + sd + 'direction of the '
            + 'photon heading toward the observer. This angle is' + sd
            + 'positive on the north side of Saturn''s rings and '
            + 'negative on the' + sd + 'south side. Its value is '
            + 'nearly constant over the duration of a' + sd
            + 'ring occultation experiment and is nearly equal to the '
            + 'ring opening' + sd + 'angle (Earth elevation angle above '
            + 'the ring plane)."')
            ,
            ('"The distance between the spacecraft (at' + sd
            + 'the SPACECRAFT EVENT TIME) and the ring-plane '
            + 'intercept point (at the' + sd + 'RING EVENT TIME)."')
            ,
            ('"Component of the displacement velocity' + sd 
            + 'of the ring-plane intercept point in the radial '
            + 'direction at the' + sd
            + 'intercept point (at the RING EVENT TIME)."')
            ,
            ('"Component of the displacement velocity' + sd
            + 'of the ring-plane intercept point in the orbital '
            + 'direction of motion' + sd + 'at the intercept point '
            + '(at the RING EVENT TIME). The orbital direction' + sd
            + 'is defined by the cross-product of a vector along '
            + 'Saturn''s pole and' + sd + 'the radial velocity at '
            + 'the intercept point, in that order."')
            ,
            ('"Fresnel scale of diffraction implied by' + sd
            + 'the occultation geometry at X-band. See Eq. (6) of '
            + 'MAROUFETAL1986.' + sd + 'The Fresnel scale at S-band is '
            + 'sqrt(11/3) times FRESNEL SCALE.' + sd + 'The Fresnel '
            + 'scale at Ka-band is sqrt(55/209) times FRESNEL SCALE."')
            ,
            ('"The radius of a sphere centered at' + sd
            + 'Saturn and is tangent to the line-of-sight from '
            + 'the spacecraft (at' + sd + 'the SPACECRAFT EVENT TIME) '
            + 'to Earth (at the OBSERVED EVENT TIME). It' + sd
            + 'identified the minimum radius of hypothetical '
            + 'spherically symmetric' + sd + 'atmosphere that is '
            + 'sensed by radio signal along its path in the' + sd
            + 'absence of refractive ray bending."')
            ,
            ('"X component of the spacecraft position' + sd
            + 'vector (at the SPACECRAFT EVENT TIME) in a '
            + 'planetocentric reference' + sd + 'frame (ux, uy, uz). '
            + 'The basis vector ux is parallel to the' + sd
            + 'projection of the line-of-sight from the spacecraft '
            + '(at the' + sd + 'SPACECRAFT EVENT TIME) to Earth '
            + '(at the OBSERVED EVENT TIME)' + sd + 'on the ring plane. '
            + 'The basis vector uz is in the' + sd
            + 'direction of Saturn''s pole. The basis vector uy '
            + 'is defined by the' + sd + 'right-hand rule."')
            ,
            ('"Y component of the spacecraft position' + sd
            + 'vector (at the SPACECRAFT EVENT TIME) in a '
            + 'planetocentric reference' + sd + 'frame (ux, uy, uz). '
            + 'See description of planetocentric reference' + sd
            + 'frame (ux, uy, uz) in the description of column '
            + 'SPACECRAFT POSITION X."')
            ,
            ('"Z component of the spacecraft position' + sd
            + 'vector (at the SPACECRAFT EVENT TIME) in a '
            + 'planetocentric reference' + sd + 'frame (ux, uy, uz). '
            + 'See description of planetocentric reference' + sd
            + 'frame (ux, uy, uz) in the description of column '
            + 'SPACECRAFT POSITION X."')
            ,
            ('"X component of the spacecraft velocity' + sd
            + 'vector (at the SPACECRAFT EVENT TIME) in a '
            + 'planetocentric reference' + sd + 'frame (ux, uy, uz). '
            + 'See description of planetocentric reference' + sd 
            + 'frame (ux, uy, uz) in the description of column '
            + 'SPACECRAFT POSITION X."')
            ,
            ('"Y component of the spacecraft velocity' + sd
            + 'vector (at the SPACECRAFT EVENT TIME) in a '
            + 'planetocentric reference' + sd + 'frame (ux, uy, uz). '
            + 'See description of planetocentric reference' + sd 
            + 'frame (ux, uy, uz) in the description of column '
            + 'SPACECRAFT POSITION X."')
            ,
            ('"Z component of the spacecraft velocity' + sd
            + 'vector (at the SPACECRAFT EVENT TIME) in a '
            + 'planetocentric reference' + sd + 'frame (ux, uy, uz). '
            + 'See description of planetocentric reference' + sd 
            + 'frame (ux, uy, uz) in the description of column '
            + 'SPACECRAFT POSITION X."')
            ,
            ('"The latitude of the position vector from' + sd
            + 'the observing DSN station to the spacecraft computed at the '
            + 'OBSERVED' + sd + 'EVENT TIME. The position vector is '
            + 'computed in the topographic reference' + sd
            + 'frame of the DSN station based on the frame kernel."')
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

def write_geo_series(rev_info, geo_inst, title, outdir, prof_dir):
    """
    This function writes a GEO series, which includes a data and label file.

    Arguments
        :rev_info (*dict*): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        :geo_inst (*class*): Instance of Geometry class
        :title (*str*): Name of the output .TAB and .LBL file, not including
                           extensions. Date in YYYYMMDD format and sequence
                           number in XXXX format will be added at the end
                           of series_name
        :outdir (*str*): Path to output directory
        :prof_dir (*str*): Direction of ring occultation for this geo_inst
    """
    outfile_tab = outdir + title.upper() + '.TAB'
    outfile_lbl = outdir + title.upper() + '.LBL'
    series_name = '"' + outfile_tab.split('/')[-1] + '"'


    # Write data file
    write_geo_series_data(geo_inst, outfile_tab)

    # Get label file information
    str_lbl = get_geo_series_info(rev_info, geo_inst, series_name, prof_dir)

    # Write label file
    pds3.pds3_write_series_lbl(str_lbl, outfile_lbl)

    return None
