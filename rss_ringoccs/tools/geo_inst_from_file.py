"""

geo_inst_from_file.py

Purpose: Create an instance with the same attributes of the Geometry class
         given a geometry file with 18 columns of same size in order:
         1) t_oet_spm_vals
         2) t_ret_spm_vals
         3) t_set_spm_vals
         4) rho_km_vals
         5) phi_rl_deg_vals
         6) phi_ora_deg_vals
         7) B_deg_vals
         8) D_km_vals
         9) rho_dot_kms_vals
         10) phi_rl_dot_kms_vals
         11) F_km_vals
         12) R_imp_km_vals
         13) rx_km_vals
         14) ry_km_vals
         15) rz_km_vals
         16) vx_kms_vals
         17) vy_kms_vals
         18) vz_kms_vals

Revisions:
    2018 Aug 02 - jfong - original
"""

import numpy as np
import pandas as pd
from write_history_dict import write_history_dict



class CreateGeoInst(object):
    """
    Purpose:
    Make a geo instance using a geo file.

    """

    def __init__(self, geo_file):
        """
        Args:
            cal_file (str): Full path name of calibration file
            rsr_inst: Instace of the RSRReader class
        """
        
        col_names = [
                't_oet_spm_vals'
                , 't_ret_spm_vals'
                , 't_set_spm_vals'
                , 'rho_km_vals'
                , 'phi_rl_deg_vals'
                , 'phi_ora_deg_vals'
                , 'B_deg_vals'
                , 'D_km_vals'
                , 'rho_dot_kms_vals'
                , 'phi_rl_dot_kms_vals'
                , 'F_km_vals'
                , 'R_imp_km_vals'
                , 'rx_km_vals'
                , 'ry_km_vals'
                , 'rz_km_vals'
                , 'vx_kms_vals'
                , 'vy_kms_vals'
                , 'vz_kms_vals'
                ]

        try:
            geo = pd.read_csv(geo_file, header=None, names=col_names)
        except FileNotFoundError:
            sys.exit('ERROR (CreateGeoInst): File not found')

        input_vars                = {"geo_file": geo_file}
        input_kwds                = {"None": None}

        self.t_oet_spm_vals       = np.asarray(geo['t_oet_spm_vals'])
        self.t_ret_spm_vals       = np.asarray(geo['t_ret_spm_vals'])
        self.t_set_spm_vals       = np.asarray(geo['t_set_spm_vals'])
        self.rho_km_vals          = np.asarray(geo['rho_km_vals'])
        self.phi_rl_deg_vals      = np.asarray(geo['phi_rl_deg_vals'])
        self.phi_ora_deg_vals     = np.asarray(geo['phi_ora_deg_vals'])
        self.B_deg_vals           = np.asarray(geo['B_deg_vals'])
        self.D_km_vals            = np.asarray(geo['D_km_vals'])
        self.rho_dot_kms_vals     = np.asarray(geo['rho_dot_kms_vals'])
        self.phi_rl_dot_kms_vals  = np.asarray(geo['phi_rl_dot_kms_vals'])
        self.F_km_vals            = np.asarray(geo['F_km_vals'])
        self.R_imp_km_vals        = np.asarray(geo['R_imp_km_vals'])
        self.rx_km_vals           = np.asarray(geo['rx_km_vals'])
        self.ry_km_vals           = np.asarray(geo['ry_km_vals'])
        self.rz_km_vals           = np.asarray(geo['rz_km_vals'])
        self.vx_kms_vals          = np.asarray(geo['vx_kms_vals'])
        self.vy_kms_vals          = np.asarray(geo['vy_kms_vals'])
        self.vz_kms_vals          = np.asarray(geo['vz_kms_vals'])
        self.kernels              = None
        self.frequency_band       = None
        self.elev_deg_vals        = None
        self.naif_toolkit_version = None
        self.beta_vals            = None
        self.B_eff_deg_vals       = None
        self.history              = write_history_dict(input_vars, input_kwds, 
                __file__)



