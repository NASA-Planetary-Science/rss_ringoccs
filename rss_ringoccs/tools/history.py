"""

:Purpose:
    Functions related to recording processing history.

:Dependencies:
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

def date_to_rev(year, doy):
    """
    Pull rev number from a table given the year and doy from a RSS
    activities file with columns for CIMS request, sequence number,
    year, doy, start earth-received time in HH:MM, end earth-received time
    in HH:MM

    Arguments
        :year (*int*): Year of occultation
        :doy (*int*): Day of year of occultation
    Keywords
        :rss_file (*str*): Table of RSS events. Default is:
            tables/RSSActivities_all_rings_only.txt

    Returns
        :rev_number (*str*): 3-digit rev number (e.g. '007')

    Note:
        #. Given default 'rss_file' location, this script must be run
            one directory from the top-level rss_ringoccs directory
    """

    rev_number_vals = np.asarray([
        "RSS_007RI_OCC003_PRIME",     "RSS_007RI_OCC004_PRIME",
        "RSS_008RI_OCC003_PRIME",     "RSS_008RI_OCC004_PRIME",
        "RSS_009RI_OCC004_PRIME",     "RSS_010RI_OCC003_PRIME",
        "RSS_010RI_OCC004_PRIME",     "RSS_011RI_OCC004_PRIME",
        "RSS_012RI_OCC003_PRIME",     "RSS_012RI_OCC004_PRIME",
        "RSS_013RI_OCC004_PRIME",     "RSS_014RI_OCC002_PRIME",
        "RSS_028RI_OCC001_PRIME",     "RSS_028RI_OCC002_PRIME",
        "RSS_044RI_OCC002_PRIME",     "RSS_046RI_OCC002_PRIME",
        "RSS_053RI_OCC002_PRIME",     "RSS_054RI_OCC002_PRIME",
        "RSS_056RI_OCC002_PRIME",     "RSS_057RI_OCC002_PRIME",
        "RSS_058RI_OCC002_PRIME",     "RSS_060RI_OCC002_PRIME",
        "RSS_063RI_OCC002_PRIME",     "RSS_064RI_OCC002_PRIME",
        "RSS_067RI_OCC002_PRIME",     "RSS_079RI_OCC002_PRIME",
        "RSS_081RI_OCC002_PRIME",     "RSS_082RI_OCC002_PRIME",
        "RSS_084RI_OCC002_PRIME",     "RSS_089RI_OCC001_PRIME",
        "RSS_123RI_OCC002_PRIME",     "RSS_123RI_OCC003_PRIME",
        "RSS_125RI_OCC002_PRIME",     "RSS_125RI_OCC003_PRIME",
        "RSS_133RI_OCC002_PRIME",     "RSS_133RI_OCC003_PRIME",
        "RSS_137RI_OCC001_PRIME",     "RSS_137RI_OCC002_PRIME",
        "RSS_167RI_OCC001_PIE",       "RSS_168RI_OCC001_PRIME",
        "RSS_169RI_OCC001_PIE",       "RSS_170RI_OCC001_PIE",
        "RSS_174RI_OCC001_PIE",       "RSS_179RI_OCC001_PRIME",
        "RSS_180RI_OCC001_PRIME",     "RSS_189RI_OCC001_PRIME",
        "RSS_190RI_OCC001_PIE",       "RSS_191RI_OCC001_PIE",
        "RSS_193RI_OCC001_PIE",       "RSS_194RI_OCC001_PRIME",
        "RSS_196RI_OCC001_PIE",       "RSS_197RI_OCC001_RSS",
        "RSS_236RI_OCC001_PIE",       "RSS_237RI_OCC001_PIE",
        "RSS_237RI_OCC002_PIE",       "RSS_238RI_OCC001_PIE",
        "RSS_247RI_OCC001_PRIME",     "RSS_248RI_OCC001_PRIME",
        "RSS_250RI_OCC001_PIE",       "RSS_251RI_OCC001_PIE",
        "RSS_253RI_OCC002_PRIME",     "RSS_255RI_OCC001_PIE",
        "RSS_256RI_OCC001_PIE",       "RSS_257RI_OCC001_PRIME",
        "RSS_266RI_OCC001_PRIME",     "RSS_268RI_OCC001_PIE",
        "RSS_270RI_OCC001_PRIME",     "RSS_273RI_PERIOCC001_PRIME",
        "RSS_273RI_INGOCC001_PRIME",  "RSS_273RI_EGROCC001_PRIME",
        "RSS_PERIOCC001_PRIME",       "RSS_274_RI_INGOCC001_PRIME",
        "RSS_274_RI_EGROCC001_PRIME", "RSS_275RI_PERIOCC001_PRIME",
        "RSS_275RI_INGOCC001_PRIME",  "RSS_275RI_EGROCC001_PRIME",
        "RSS_276RI_OCC001_PRIME",     "RSS_278RI_PERIOCC001_PRIME",
        "RSS_278RI_CRDOCC001_PRIME",  "RSS_280RI_PERIOCC001_PRIME",
        "RSS_280RI_CRDOCC001_PRIME",  "RSS_282RI_OCC001_PRIME",
        "RSS_284RI_PERIOCC001_PRIME", "RSS_284RI_CRDOCC001_PRIME"
    ])
    year_vals = np.asarray([
        2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005, 2005,
        2005, 2006, 2006, 2007, 2007, 2007, 2007, 2008, 2008, 2008, 2008,
        2008, 2008, 2008, 2008, 2008, 2008, 2008, 2008, 2009, 2009, 2010,
        2010, 2010, 2010, 2010, 2010, 2012, 2012, 2012, 2012, 2012, 2013,
        2013, 2013, 2013, 2013, 2013, 2013, 2013, 2013, 2016, 2016, 2016,
        2016, 2016, 2016, 2016, 2016, 2016, 2017, 2017, 2017, 2017, 2017,
        2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017, 2017,
        2017, 2017, 2017, 2017, 2017, 2017, 2017
    ])
    doy_vals = np.asarray([
        123, 123, 141, 141, 159, 177, 177, 196, 214, 214, 232, 248, 258,
        259, 130, 162, 337, 353,  15,  27,  39,  62,  92, 102, 130, 217,
        232, 239, 254, 291, 359, 360,  26,  27, 169, 170, 245, 245, 156,
        180, 204, 225, 315,  18,  31, 130, 140, 151, 175, 187, 220, 244,
        158, 182, 182, 206, 307, 317, 333, 340, 354,   3,  10,  17,  81,
         96, 110, 129, 129, 129, 135, 135, 136, 142, 142, 142, 148, 161,
        161, 174, 174, 187, 200, 200
    ])

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

def rev_to_occ_info(rev):
    """
    Pull occultation direction from a text file given rev.

    Arguments
        :rev (*str*):  Revolution/orbit number in 'XXX' format

    Keyword Arguments
        :sroc_info_file (*str*): Path to csv file with columns: rev number,
            occultation direction, planetary occultation flag. Default values
            come from tables/list_of_sroc_dir_all_events.txt.

    Returns
        :occ_dir (*str*): Occultation direction (over entire, I&E, occultation)
                       This is not to be confused with profile direction.

    Note:
        #. Given default 'sroc_info_file' location, this script must be run
            one directory from the top-level rss_ringoccs directory
    """

    rev_str_list = [
        "007", "008", "009", "010", "011", "012", "013", "014", "028",
        "044", "046", "053", "054", "056", "057", "058", "060", "063",
        "064", "067", "079", "081", "082", "084", "089", "123", "125",
        "133", "137", "167", "168", "169", "170", "174", "179", "180",
        "189", "190", "191", "193", "194", "196", "197", "236", "237",
        "238", "247", "248", "250", "251", "253", "255", "256", "257",
        "266", "268", "270", "273", "274", "275", "276", "278", "280",
        "282", "284"
    ]
    occ_dir_list = [
        "BOTH",    "BOTH",    "EGRESS",  "BOTH",    "EGRESS",  "BOTH",
        "EGRESS",  "INGRESS", "INGRESS", "EGRESS",  "INGRESS", "BOTH",
        "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",
        "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",
        "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",
        "INGRESS", "INGRESS", "BOTH",    "INGRESS", "INGRESS", "BOTH",
        "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH",
        "BOTH",    "BOTH",    "BOTH",    "EGRESS",  "BOTH",    "BOTH",
        "BOTH",    "BOTH",    "BOTH",    "EGRESS",  "EGRESS",  "EGRESS",
        "INGRESS", "INGRESS", "INGRESS", "BOTH",    "BOTH",    "BOTH",
        "BOTH",    "BOTH",    "BOTH",    "BOTH",    "BOTH"
    ]

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
    rssocc_version = '1.3-beta'

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
