"""
################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or              #
#   modify it under the terms of the GNU General Public License as published   #
#   by the Free Software Foundation, either version 3 of the License, or       #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful             #
#   but WITHOUT ANY WARRANTY# without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Purpose:                                                                   #
#       Creates an LBL file from Geo data.                                     #
################################################################################
#   Author: Jolene Fong                                                        #
#   Date:   2018/07/23                                                         #
################################################################################
"""
# PDS names are all caps. Pylint doesn't like this.
# pylint: disable = invalid-name
import time
import numpy
from . import pds3_write_series as pds3

def write_geo_series_data(geo_inst, fmt, out_file):
    """
    This writes a GEO data file.

    Args:
        geo_inst (class): Instance of Geometry class
        fmt (str): Format string
        out_file (str): Output file name, including path.
    """
    fmt_comma = fmt+','
    format_str = fmt_comma * 17 + fmt + '%s'
    npts = len(geo_inst.t_oet_spm_vals)

    print('\nWriting GEO data to: ', out_file, '\n')

    with open(out_file, 'w', encoding = "utf8") as f:

        for n in range(npts):
            f.write(
                format_str % (
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
                '\r\n'
            )
        )

    f.close()

def get_geo_series_info(rev_info, geo_inst, series_name, prof_dir):
    """
    This returns the information needed to write a GEO label file.

    Args:
        rev_info (dict): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        geo_inst (class): Instance of Geometry class
        series_name (str): Name of the output .TAB and .LBL file, not including
                           extensions. Date in YYYYMMDD format will be added
                           onto series_name
        prof_dir (str): Direction of ring occultation for this geo_inst

    Outputs:
        str_lbl (dict): Dictionary with keys: string_delimiter,
                        alignment_column, series_alignment_column,
                        keywords_value, keywords_NAIF_TOOLKIT_VERSION,
                        description, keywords_series, object_keys,
                        object_values, history

    Notes:
        [1] This is a reproduction of GEO label files within
            Cassini_RSS_Ring_Profiles_2018_Archive, with minor edits.
        [2] The format of each data entry is hardcoded within "nchar".
    """
    # Get current time in ISOD format
    current_time_ISOD = time.strftime("%Y-%j") + 'T' + time.strftime("%H:%M:%S")

    # Values for number of columns, number of bytes per column,
    #   number of bytes per column delimiter, number of bytes allocated to
    #   special characters per column
    ncol = 18
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
    start_bytes_list = numpy.cumsum(new_bytes_list)
    start_bytes_list = start_bytes_list[:-1]

    col_num = list(range(1, ncol+1))

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
    geo_kernels = ['"'+x.split('/')[-1]+'"' for x in geo_kernels]

    PDS_VERSION_ID = 'PDS3'
    RECORD_TYPE = 'FIXED_LENGTH'
    RECORD_BYTES = pds3.get_record_bytes(ncol, nchar, ndelim, nspecial)
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

    sd = '|'

    # Description of file, with sd marking end of lines
    FILE_DESCRIPTION = ('"This file contains parameters' + sd
            + 'pertinent to the ring occultation observation'
            + 'geometry needed for' + sd + 'reconstruction of'
            + 'diffraction effects and for other data analysis'
            + sd + 'purposes. Light-time effects have been'
            + 'accounted for in computing' + sd + 'these parameters. '
            + 'For more details, see MAROUFETAL1986 and Appendix' + sd
            + 'B of ROSEN1989.' + sd
            + ' ' + sd
            + 'Please address any inquiries to:' + sd
            + 'Richard G. French' + sd
            + 'Astronomy Department, Wellesley College' + sd
            + 'Wellesley, MA 02481-8203' + sd
            + '(781) 283-3747' + sd
            + 'rfrench@wellesley.edu"')

    # Extract history information
    HIST_USER_NAME = geo_inst.history['User Name']
    HIST_HOST_NAME = geo_inst.history['Host Name']
    HIST_RUN_DATE = geo_inst.history['Run Date']
    HIST_PYTHON_VERSION = geo_inst.history['Python Version']
    HIST_OPERATING_SYSTEM = geo_inst.history['Operating System']
    HIST_SOURCE_DIR = geo_inst.history['Source Directory']
    HIST_SOURCE_FILE = geo_inst.history['Source File']
    HIST_INPUT_VARIABLES = geo_inst.history['Input Variables']
    HIST_INPUT_KEYWORDS = geo_inst.history['Input Keywords']
    HIST_description = ('This is a detailed record of the' + sd
                    + 'processing steps used to generate this file.')

    HISTORY_dict = {
            'key_order': ['User Name', 'Host Name', 'Run Date',
                        'Python Version', 'Operating System',
                        'Source Directory','Source File',
                        'Input Variables', 'Input Keywords']
            , 'hist name': 'Geometry history'
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


    # Write series information
    SERIES = 'SERIES'
    SERIES_NAME = '"OCCULTATION GEOMETRY"'
    SERIES_INTERCHANGE_FORMAT = 'ASCII'
    SERIES_COLUMNS = '18'
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
            ]

    n_objects = len(object_names)
    data_types = ['ASCII_REAL'] * n_objects
    # NOTE: this is hardcoded in!! and a bad way to make a list of reps
    formats = ['"F32.16"'] * 18
    units = ['"SECOND"', '"SECOND"', '"SECOND"',
            '"KILOMETER"', '"DEGREE"', '"DEGREE"', '"DEGREE"',
            '"KILOMETER"', '"KM/SEC"', '"KM/SEC"',
            '"KILOMETER"', '"KILOMETER"', '"KILOMETER"', '"KILOMETER"',
            '"KILOMETER"', '"KM/SEC"', '"KM/SEC"', '"KM/SEC"']

    # only variables with TIME in name have reference times
    # a null string indicates that this keyword is absent for this object
    # NOTE: this is a horrible way of creating a list! (same as units)
    es = ''
    reference_times = [OBJECT_REFERENCE_TIME, OBJECT_REFERENCE_TIME,
            OBJECT_REFERENCE_TIME, es, es, es, es, es, es, es, es, es,
            es, es, es, es, es, es]

    object_descriptions = [
            ('"The instant at which photons were' + sd
            + 'received at the DSN receiving station, given in '
            + 'elapsed seconds' + sd + 'after the moment specified '
            + 'by REFERENCE_TIME."')
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
            + 'the moment specified by' + sd + 'REFERENCE_TIME."')
            ,
            ('"Radial distance from the center of' + sd
            + 'Saturn to the ring-plane intercept point at '
            + 'the RING EVENT TIME."')
            ,
            ('"Inertial longitude on the ring plane' + sd
            + 'of the ring-plane intercept point at the RING EVENT TIME."')
            ,
            ('"Measured at the ring-plane intercept' + sd
            + 'point, starting from the direction of a photon heading to '
            + 'the' + sd + 'observer, and ending at the direction '
            + 'of a local radial vector.' + sd + 'This angle is '
            + 'projected into the ring plane and measured in the' + sd
            + 'prograde direction. As seen from the observer, it equals 90 '
            + sd + 'degrees along the right ansa and 270 degrees '
            + 'along the left ansa.' + sd + 'Values range from 0 '
            + 'to 360 degrees. This convention for observed' + sd
            + 'ring azimuth differs from that adopted in '
            + 'MAROUFETAL1986 by 180' + sd + 'degrees."')
            ,
            ('"The angle measured at the ring' + sd
            + 'intercept point, starting from the ring plane '
            + 'and ending in the' + sd + 'direction of the '
            + 'photon heading toward the observer. This angle is' + sd
            + 'positive on the north side of Saturn''s rings and '
            + 'negative on the' + sd + 'south side."')
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
            + 'MAROUFETAL1986.' + sd + 'The Fresnel scale at S-band is'
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

    Args:
        rev_info (dict): Dictionary with keys: rsr_file, band, year, doy, dsn
                         occ_dir, planetary_occ_flag, rev_num
        geo_inst (class): Instance of NormDiff class
        title (str): Name of the output .TAB and .LBL file, not including
                           extensions. Date in YYYYMMDD format will be added
                           onto series_name
        outdir (str): Path to output directory
        prof_dir (str): Direction of ring occultation for this geo_inst

    Notes:
        [1] Data entry format of %32.16F is hardcoded.
        [2] A data and label file will be output into the input "outdir"
            directory, with filenames, *YYYYMMDD.TAB and *YYYYMMDD.LBL,
            respectively, where * is "title".
    """
    current_time = time.strftime("_%Y%m%d") #-%H%M%S")
    outfile_tab = outdir + title.upper() + current_time + '.TAB'
    outfile_lbl = outdir + title.upper() + current_time + '.LBL'
    series_name = '"' + outfile_tab.split('/')[-1] + '"'

    fmt = '%32.16F' # format for float entry


    # Write data file
    write_geo_series_data(geo_inst, fmt, outfile_tab)

    # Get label file information
    str_lbl = get_geo_series_info(rev_info, geo_inst, series_name, prof_dir)

    # Write label file
    pds3.pds3_write_series_lbl(str_lbl, outfile_lbl)
