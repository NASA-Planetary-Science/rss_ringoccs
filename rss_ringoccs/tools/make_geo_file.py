"""

jwf_make_geo_file.py

"""
import numpy as np

def make_geo_file(geo_inst, geo_file):
	format_str = ['%20.9E'] * 18
	input_text = np.c_[geo_inst.t_oet_spm_vals,
					geo_inst.t_ret_spm_vals,
					geo_inst.t_set_spm_vals,
					geo_inst.rho_km_vals,
					geo_inst.phi_rl_deg_vals,
					geo_inst.phi_ora_deg_vals,
					geo_inst.B_deg_vals,
					geo_inst.D_km_vals,
					geo_inst.rho_dot_kms_vals,
					geo_inst.phi_rl_dot_kms_vals,
					geo_inst.F_km_vals,
					geo_inst.R_imp_km_vals,
					geo_inst.rx_km_vals,
					geo_inst.ry_km_vals,
					geo_inst.rz_km_vals,
					geo_inst.vx_kms_vals,
					geo_inst.vy_kms_vals,
					geo_inst.vz_kms_vals]
	np.savetxt(geo_file, input_text, delimiter = ',', fmt=format_str)
	return None

