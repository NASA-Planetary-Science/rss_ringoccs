"""

get_rev_info.py

Purpose: Create dictionary of rev information needed for LBL files.

Revisions:
    2018 Jul 18 - jfong - original
                        - move rev_to_occ_info.py here


"""

import pandas as pd
import sys
import pdb

def get_rev_info(rsr_inst, rev):
    """
    This returns a dictionary with information related to the ring occultation.

    Args:
        rev (str): Revolution/orbit number in 'XXX' format

    Output:
        rev_info (dict): Dictionary with keys: rsr_file, band, year, doy
                                                dsn, rev, occ_dir, 
                                                planetary_occ_flag
    Notes:
        [1] This is meant to be run from the pipeline directory because of the
            relative path to the table file.
        [2] This will only work for occultations BEFORE USO failure (up to 137)

    """

    occ_dir, planetary_occ_flag = rev_to_occ_info(rev)

    rev_info = {
            "rsr_file":   rsr_inst.rsr_file.split('/')[-1]
            , "band":     '"'+str(rsr_inst.band).split("'")[1]+'"'
            , "year":     str(rsr_inst.year)
            , "doy":      str(rsr_inst.doy)
            , "dsn":      str(rsr_inst.dsn)
            , "occ_dir":  occ_dir
            , "planetary_occ_flag": planetary_occ_flag
            , "rev_num": rev
            }

    return rev_info



def rev_to_occ_info(rev):
    """
    Pull occultation direction from a text file given rev.

    Args:
        rev (str):  Revolution/orbit number in 'XXX' format

    Output:
        occ_dir (str): Occultation direction (over entire, I&E, occultation)
                       This is not to be confused with profile direction.
        planetary_occ_flag (str): Flag for whether planet (Saturn) was
                                  occulted during event.
    Notes:
        [1] This is meant to be run from the pipeline directory because of the
            relative path to the table file.
        [2] This will only work for occultations BEFORE USO failure (up to 137)


    """
    csv_file = '../tables/list_of_sroc_dir_before_USOfailure.txt'

    occ_dir_table = pd.read_csv(csv_file, header=None, skiprows=1, dtype=str)
    rev_str_list = list(occ_dir_table[0])
    occ_dir_list = list(occ_dir_table[1])
    planet_flag_list = list(occ_dir_table[2])

    try:
        ind = rev_str_list.index(rev)
    except ValueError as err:
        sys.exit('WARNING (rev_to_occ_info): rev '+ str(err))

    occ_dir = '"' + occ_dir_list[ind] + '"'
    planetary_occ_flag = '"' + planet_flag_list[ind] + '"'

    return occ_dir, planetary_occ_flag



