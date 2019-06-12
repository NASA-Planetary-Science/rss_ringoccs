"""
history.py

Purpose:
    Functions related to recording processing history.

Dependencies:
    #. sys
    #. time
    #. os
    #. platform
    #. pandas
    #. numpy
"""
import sys
import time
import platform
import os
import numpy as np
import pandas as pd

def date_to_rev(year, doy,
        rss_file='../tables/RSSActivities_all_rings_only.txt'):
    """
    Purpose:
        Pull rev number from a table given the year and doy from a RSS
        activities file with columns for CIMS request, sequence number,
        year, doy, start earth-received time in HH:MM, end earth-received time
        in HH:MM

    Arguments:
        :year (*int*): Year of occultation
        :doy (*int*): Day of year of occultation

    Returns:
        :rev_number (*str*): 3-digit rev number (e.g. '007')

    Note:
        #. Given default 'rss_file' location, this script must be run
            one directory from the top-level rss_ringoccs directory
    """

    date_to_rev_table = pd.read_csv(rss_file, header=None, skiprows=1)
    rev_number_vals = np.asarray(date_to_rev_table[0])
    year_vals = np.asarray(date_to_rev_table[2])
    doy_vals = np.asarray(date_to_rev_table[3])
    year_where = (year_vals == year)
    doy_where = (doy_vals == doy)
    if (doy_where == False).all():
        possible_doy = list(doy_vals[year_where])
        ind = possible_doy.index(min(possible_doy, key=lambda x:abs(x-doy)))
        doy_where = (doy_vals == possible_doy[ind])
    year_and_doy_where = year_where & doy_where
    rev_number = rev_number_vals[year_and_doy_where][0][4:7]

    return rev_number

def get_rev_info(rsr_inst):
    """
    This returns a dictionary with information related to the ring occultation.

    Arguments
        :rsr_inst (*class*): Instance of RSRReader class

    Returns:
        :rev_info (*dict*): Dictionary with keys: rsr_file, band, year, doy
                            dsn, rev, occ_dir, planetary_occ_flag
    """

    rev = date_to_rev(rsr_inst.year, rsr_inst.doy)

    occ_dir = rev_to_occ_info(rev)

    rev_info = {
            "rsr_file":   rsr_inst.rsr_file.split('/')[-1]
            , "band":     '"'+str(rsr_inst.band)+'"'
            , "year":     str(rsr_inst.year)
            , "doy":      str(rsr_inst.doy).zfill(3)
            , "dsn":      str(rsr_inst.dsn)
            , "occ_dir":  occ_dir
            , "rev_num": rev
            }
           # , "planetary_occ_flag": planetary_occ_flag

    return rev_info

def rev_to_occ_info(rev,
        sroc_info_file='../tables/list_of_sroc_dir_all_events.txt'):
    """
    Pull occultation direction from a text file given rev.

    Arguments
        :rev (*str*):  Revolution/orbit number in 'XXX' format

    Keyword Arguments
        :sroc_info_file (*str*): Path to csv file with columns: rev number,
            occultation direction, planetary occultation flag

    Returns
        :occ_dir (*str*): Occultation direction (over entire, I&E, occultation)
                       This is not to be confused with profile direction.
        :planetary_occ_flag (*str*): Flag for whether planet (Saturn) was
                                  occulted during event.
    Note:
        #. Given default 'sroc_info_file' location, this script must be run
            one directory from the top-level rss_ringoccs directory
    """

    occ_dir_table = pd.read_csv(sroc_info_file, header=None, skiprows=1,
            dtype=str)
    rev_str_list = list(occ_dir_table[0])
    occ_dir_list = list(occ_dir_table[1])

    try:
        ind = rev_str_list.index(rev)
    except ValueError:
        print('(rev_to_occ_info()): Rev not found in sroc_info_file!')
        return '"BOTH"'#,'"Y"'

    occ_dir = '"' + occ_dir_list[ind] + '"'

    return occ_dir

def write_history_dict(input_vars, input_kwds, source_file, add_info=None):
    """
    This creates a dictionary of processing history for an instance.

    Arguments:
        :input_vars (*dict*): Dictionary of all input variables to the
                                instance.
        :input_kwds (*dict*): Dictionary of all input keywords to the instance.
        :source_file (*str*):  Full path to the script used to run the instance.

    Keyword Arguments:
        :add_info (*dict*): Dictionary of additional info

    Returns:
        :history (*dict*): Dictionary with keys: "User Name", "Host Name",
                            "Run Date", "Python Version", "Operating System",
                            "Source File", "Positional Args",
                            "Keyword Args"

    """

    user_name = os.getlogin()
    host_name = os.uname()[1]
    run_date = time.ctime() + ' ' + time.tzname[0]
    python_version = platform.python_version()
    operating_system = os.uname()[0]
    src_dir = source_file.rsplit('/',1)[0] +'/'
    src_file = source_file.split('/')[-1]
    rssocc_version = '1.1'

    history = {
            "rss_ringoccs Version": rssocc_version,
            "User Name": user_name,
            "Host Name": host_name,
            "Run Date": run_date,
            "Python Version": python_version,
            "Operating System": operating_system,
            "Source Directory": src_dir,
            "Source File": src_file,
            "Positional Args": input_vars,
            "Keyword Args": input_kwds
            }

    if add_info:
        history["Additional Info"] = add_info
    else:
        history["Additional Info"] = ''
    return history
