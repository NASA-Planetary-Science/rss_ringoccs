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
#       Pull rev number from a table, given year and doy.                      #
################################################################################
#   Author: Jolene Fong                                                        #
#   Date:   2018/07/12                                                         #
################################################################################
"""
# Pylint claims year_where = (year_vals == year) uses superfluous parenthesis.
# But year_where = year_vals == year is confusing syntax. Ignore this warning.
# pylint: disable = superfluous-parens
import numpy
import pandas

def date_to_rev(year, doy):
    """
        Function:
            date_to_rev
        Purpose:
            Pull rev number from a table given the year and doy. There are
            often multiple rows matching the year and doy, but these
            multiple rows always have the same rev number, so that's okay.
            Just take the first one.
        Args:
            year (int):
                Year of occultation
            doy (int):
                Day of year of occultation
        Output:
            rev_number:
                The Rev number from the given date.
        Notes:
            Must be run within pipeline directory because
            of the relative path of the data file.
    """

    date_to_rev_table = pandas.read_csv(
        '../tables/RSSActivities_before_USOfailure_rings_only.txt',
        header=None, skiprows=1
    )

    rev_number_vals = numpy.asarray(date_to_rev_table[0])
    year_vals = numpy.asarray(date_to_rev_table[2])
    doy_vals = numpy.asarray(date_to_rev_table[3])
    year_where = (year_vals == year)
    doy_where = (doy_vals == doy)
    year_and_doy_where = year_where & doy_where
    rev_number = rev_number_vals[year_and_doy_where][0][4:7]

    return rev_number
