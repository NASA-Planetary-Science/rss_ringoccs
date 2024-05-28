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
#       Create an instance with the same attributes as the Geometry class      #
#       given a geometry file with 18 columns of the same size in order.       #
################################################################################
#                                 DEPENDENCIES                                 #
################################################################################
#   1.) numpy:                                                                 #
#           Used for creating arrays from the data.                            #
#   2.) pandas:                                                                #
#           Extracts the data from the CSV file.                               #
################################################################################
#   Author: Jolene Fong                                                        #
#   Date:   2018/08/02                                                         #
################################################################################
"""

# These geometry files have a lot of variables. Pylint may complain about this.
# pylint: disable = too-many-instance-attributes

# The naming convention also uses upper case for certain variables like the
# ring opening angle and Fresnel scale. Ignore these warnings too.
# pylint: disable = invalid-name
import numpy
import pandas
from .write_history_dict import write_history_dict

class GeoInstanceFromFile:
    """
        Class:
            GeoInstanceFromFile
        Purpose:
            Make a geo instance using a geo file.
        Arguments:
            geo_file (str):
                A geometry CSV file. It must have the following data, in order:
                    t_oet_spm_vals
                    t_ret_spm_vals
                    t_set_spm_vals
                    rho_km_vals
                    phi_rl_deg_vals
                    phi_ora_deg_vals
                    B_deg_vals
                    D_km_vals
                    rho_dot_kms_vals
                    phi_rl_dot_kms_vals
                    F_km_vals
                    R_imp_km_vals
                    rx_km_vals
                    ry_km_vals
                    rz_km_vals
                    vx_kms_vals
                    vy_kms_vals
                    vz_kms_vals
        History:
            2018/08/02: Jolene Fong
                Original draft.
            2024/05/28: Ryan Maguire
                Added license, cleaned up, and added to graveyard directory.
    """
    def __init__(self, geo_file):
        """
            Args:
                geo_file (str):
                    Full path name of the geometry file.
        """

        col_names = [
            "t_oet_spm_vals",
            "t_ret_spm_vals",
            "t_set_spm_vals",
            "rho_km_vals",
            "phi_rl_deg_vals",
            "phi_ora_deg_vals",
            "B_deg_vals",
            "D_km_vals",
            "rho_dot_kms_vals",
            "phi_rl_dot_kms_vals",
            "F_km_vals",
            "R_imp_km_vals",
            "rx_km_vals",
            "ry_km_vals",
            "rz_km_vals",
            "vx_kms_vals",
            "vy_kms_vals",
            "vz_kms_vals"
        ]

        self.set_attributes_to_none()

        geo = pandas.read_csv(geo_file, header = None, names = col_names)

        self.t_oet_spm_vals = numpy.asarray(geo["t_oet_spm_vals"])
        self.t_ret_spm_vals = numpy.asarray(geo["t_ret_spm_vals"])
        self.t_set_spm_vals = numpy.asarray(geo["t_set_spm_vals"])
        self.rho_km_vals = numpy.asarray(geo["rho_km_vals"])
        self.phi_rl_deg_vals = numpy.asarray(geo["phi_rl_deg_vals"])
        self.phi_ora_deg_vals = numpy.asarray(geo["phi_ora_deg_vals"])
        self.B_deg_vals = numpy.asarray(geo["B_deg_vals"])
        self.D_km_vals = numpy.asarray(geo["D_km_vals"])
        self.rho_dot_kms_vals = numpy.asarray(geo["rho_dot_kms_vals"])
        self.phi_rl_dot_kms_vals  = numpy.asarray(geo["phi_rl_dot_kms_vals"])
        self.F_km_vals = numpy.asarray(geo["F_km_vals"])
        self.R_imp_km_vals = numpy.asarray(geo["R_imp_km_vals"])
        self.rx_km_vals = numpy.asarray(geo["rx_km_vals"])
        self.ry_km_vals = numpy.asarray(geo["ry_km_vals"])
        self.rz_km_vals = numpy.asarray(geo["rz_km_vals"])
        self.vx_kms_vals = numpy.asarray(geo["vx_kms_vals"])
        self.vy_kms_vals = numpy.asarray(geo["vy_kms_vals"])
        self.vz_kms_vals = numpy.asarray(geo["vz_kms_vals"])

        self.history = self.get_history(geo_file)

    def set_attributes_to_none(self):
        """
            Sets all attributes to None.
        """
        self.t_oet_spm_vals = None
        self.t_ret_spm_vals = None
        self.t_set_spm_vals = None
        self.rho_km_vals = None
        self.phi_rl_deg_vals = None
        self.phi_ora_deg_vals = None
        self.B_deg_vals = None
        self.D_km_vals = None
        self.rho_dot_kms_vals = None
        self.phi_rl_dot_kms_vals = None
        self.F_km_vals = None
        self.R_imp_km_vals = None
        self.rx_km_vals = None
        self.ry_km_vals = None
        self.rz_km_vals = None
        self.vx_kms_vals = None
        self.vy_kms_vals = None
        self.vz_kms_vals = None
        self.kernels = None
        self.frequency_band = None
        self.elev_deg_vals = None
        self.naif_toolkit_version = None
        self.beta_vals = None
        self.B_eff_deg_vals = None

    def get_history(self, geo_file):
        """
            Returns the history dictionary for the given geo file.
        """
        input_vars = {"geo_file": geo_file}
        input_kwds = {"None": None}
        return write_history_dict(input_vars, input_kwds, __file__)
