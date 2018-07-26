"""

date_to_rev.py

Purpose: Pull rev number from a table, given year and doy.

Revisions:
    2018 Jul 12 - jfong - copied from gsteranka's RSS_loop_v1
                        - removed 'Rev' from rev_number
"""
import numpy as np 
import pandas as pd


def date_to_rev(year, doy):
    """
    Pull rev number from a table given the year and doy. There are often
    multiple rows matching the year and doy, but these multiple rows always
    have the same rev number, so that's okay. Just take the first one.

    Args:
        year (int): Year of occultation
        doy (int): Day of year of occultation

    Note:
        [1] must be run within pipeline directory because of the relative path
            of the data file
    """

    date_to_rev_table = pd.read_csv(
        '../tables/RSSActivities_before_USOfailure_rings_only.txt',
        header=None, skiprows=1)
    rev_number_vals = np.asarray(date_to_rev_table[0])
    year_vals = np.asarray(date_to_rev_table[2])
    doy_vals = np.asarray(date_to_rev_table[3])
    year_where = (year_vals == year)
    doy_where = (doy_vals == doy)
    year_and_doy_where = year_where & doy_where
    rev_number = rev_number_vals[year_and_doy_where][0][4:7]

    return rev_number
