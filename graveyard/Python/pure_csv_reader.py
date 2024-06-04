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
#       Extracts variables from a CSV file.                                    #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018                                                               #
################################################################################
"""
# Disable a few of the annoying pylint warnings. The variable names use
# capital letters, and the CSV files have more than 7 variables associated to
# them.
# pylint: disable = invalid-name
# pylint: disable = too-many-instance-attributes
import numpy
import pandas

class PureCSVReader:
    """
        Class:
            PureCSVReader
        Purpose:
            Reads data from a CSV file and does not alter them in any way.
            That is, no interpolation is done, no smoothing, etc.
        Attribute:
            rho_km_vals (numpy.ndarray):
                The radii of the points in the data set.
            phi_rad_vals (numpy.ndarray):
                The azimuthal angles for the points.
            p_norm_vals (numpy.ndarray):
                The normalized power as a function of rho.
            phase_rad_vals (numpy.ndarray):
                The complex angle, or phase, for the points.
            f_sky_hz_vals (numpy.ndarray):
                The frequency of the data as a function of rho.
            B_rad_vals (numpy.ndarray):
                The ring opening angle as a function of rho.
            D_km_vals (numpy.ndarray):
                The spacecraft-to-ring distance.
            rho_dot_kns_vals (numpy.ndarray):
                The time derivative of rho as a function of rho.
    """
    def __init__(self, dat):
        """
            Method:
                __init__:
            Purpose:
                Constructs the class.
            Arguments:
                dat (str):
                    Path to the CSV file.
        """
        if not isinstance(dat, str):
            raise TypeError("Text file must be a string.")

        data_file = pandas.read_csv(dat)
        self.set_variables_from_data_file(data_file)
        self.check_variable_lengths()

    def set_variables_from_data_file(self, data_file):
        """
            Method:
                set_variables_from_data_file
            Purpose:
                Creates numpy arrays from the data in a CSV file.
            Arguments:
                data_file:
                    A pandas object holding the CSV data.
        """
        self.rho_km_vals = numpy.array(data_file.rho_km_vals)
        self.phase_rad_vals = numpy.array(data_file.phi_rad_vals)
        self.p_norm_vals = numpy.array(data_file.p_norm_vals)
        self.phi_rad_vals = numpy.array(data_file.phase_rad_vals)
        self.f_sky_hz_vals = numpy.array(data_file.f_sky_hz_vals)
        self.B_rad_vals = numpy.array(data_file.B_rad_vals)
        self.D_km_vals = numpy.array(data_file.D_km_vals)
        self.rho_dot_kms_vals = numpy.array(data_file.rho_dot_kms_vals)

    def check_variable_lengths(self):
        """
            Method:
                check_variable_lengths
            Purpose:
                Checks the sizes of the numpy arrays.
        """

        # All of the variables to check.
        variable_list = [
            self.phi_rad_vals,
            self.p_norm_vals,
            self.phase_rad_vals,
            self.B_rad_vals,
            self.f_sky_hz_vals,
            self.D_km_vals,
            self.rho_dot_kms_vals
        ]

        # The lengths should match the length of rho.
        array_length = numpy.size(self.rho_km_vals)

        # Loop through and check the data.
        for variable in variable_list:
            if numpy.size(variable) != array_length:
                raise IndexError("Array lengths are not equal.")
