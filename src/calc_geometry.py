
'''

jwf_calc_geometry.py
Purpose: Class Geom whose attributes are the occultation geometry parameters.
Revisions:
	2018 Feb 12 - jfong - original
	2018 Feb 17 - jfong - compatible with jwf_read_rsr_v2
	2018 Feb 22 - jfong - add test_file kwd that compares results to given file
	2018 Mar 02 - jfong - test_file kwd plots fractional error

'''

import spiceypy as spice
import numpy as np
from jwf_rsr_reader_v2 import RSRReader
import matplotlib.pyplot as plt
import time


# Global variables -- will probably be in some main program?
SPACECRAFT              = 'Cassini'
PLANET                  = 'Saturn'
SATURN_NAIF_ID          = 699
SATURN_NAIF_RADII       = 'BODY699_RADII'

class Geom(object):

	"""This is an object that calculates occultation geometry from event time.
	"""

	def __init__(self, rsr_header, kernels, pt_per_sec=None, test_file=None, 
			write_file=None):
		"""This calculates occultation geometry as attributes.

		Args:
			rsr_header(obj): 
				Header of RSR file containing SPM and F_sky_hz
			kernels(list(str)): 
				List of SPICE kernels that cover event time.
			pt_per_sec(int, optional): 
				Number of points calculated per second for all attributes.
				This will default to 1.
			test_file(str, optional): 
				Filepath to a GEO data file to compare results to
			write_file(str, optional): 
				Filepath of output GEO data file
		"""

		self.t_OET_spm_vals			= None
		self.RET_ET_vals			= None
		self.SET_ET_vals			= None
		self.rho_km_vals			= None
		self.phi_RL_deg_vals		= None
		self.phi_ORA_deg_vals		= None
		self.B_deg_vals				= None
		self.D_km_vals				= None
		self.rho_dot_kms_vals		= None
		self.phi_RL_dot_kms_vals	= None
		self.F_km_vals				= None
		self.R_imp_km_vals			= None
		self.rx_km_vals				= None
		self.ry_km_vals				= None
		self.rz_km_vals				= None
		self.vx_kms_vals			= None
		self.vy_kms_vals			= None
		self.vz_kms_vals			= None

		if pt_per_sec == None:
			pt_per_sec = 1.

		
		# Create SPM array with defined points per seconc
		spm_start = rsr_header.SPM_vals[0]
		spm_end = rsr_header.SPM_vals[-1]

		step = 1./pt_per_sec
		self.t_OET_spm_vals = np.arange(spm_start, spm_end, step)
		
		F_sky_hz_vals = np.interp(self.t_OET_spm_vals, rsr_header.f_SPM, 
				rsr_header.f_sky_pred)


		spice.kclear()
		for kernel in kernels:
			spice.furnsh(kernel)

		year	= rsr_header.year
		DOY		= rsr_header.DOY
		DSN		= 'DSS-'+str(rsr_header.DSN)
		band	= rsr_header.band

		rho_km_vals, spoint_J2000, R_sc2dsn_km_vals = self.calc_dynam_frame( 
				self.t_OET_spm_vals, year, DOY, DSN)


		phi_RL_deg_vals, phi_ORA_deg_vals = self.calc_phi(spoint_J2000, 
				R_sc2dsn_km_vals)

		D_km_vals = self.calc_D(self.RET_ET_vals, self.SET_ET_vals)

		B_deg_vals = self.calc_B(self.SET_ET_vals, DSN)

		F_km_vals = self.calc_F(D_km_vals, F_sky_hz_vals, B_deg_vals, 
				phi_ORA_deg_vals)

		R_sc_km_vals, R_sc_dot_kms_vals, R_imp_km_vals = \
				self.calc_SC_state(self.SET_ET_vals, DSN)

		t_RET_spm_vals = self.et2spm(self.RET_ET_vals)
		t_SET_spm_vals = self.et2spm(self.SET_ET_vals)

		rho_dot_kms_vals = rho_km_vals - np.roll(rho_km_vals, 1)
		phi_RL_rad_vals = [x*spice.rpd() for x in phi_RL_deg_vals]
		phi_RL_dot_rads_vals = phi_RL_rad_vals - np.roll(phi_RL_rad_vals, 1)
		phi_RL_dot_kms_vals	= [a*b for a,b in zip(rho_km_vals, 
						phi_RL_dot_rads_vals)]

		self.t_RET_SPM_vals = t_RET_spm_vals
		self.t_SET_SPM_vals = t_SET_spm_vals
		self.rho_km_vals = rho_km_vals
		self.phi_RL_deg_vals = phi_RL_deg_vals
		self.phi_ORA_deg_vals = phi_ORA_deg_vals
		self.D_km_vals = D_km_vals
		self.B_deg_vals = B_deg_vals
		self.rho_dot_kms_vals = rho_dot_kms_vals
		self.phi_RL_dot_kms_vals = phi_RL_dot_kms_vals
		self.F_km_vals = F_km_vals
		self.R_imp_km_vals = R_imp_km_vals
		self.rx_km_vals = np.stack(R_sc_km_vals)[:,0]
		self.ry_km_vals = np.stack(R_sc_km_vals)[:,1]
		self.rz_km_vals = np.stack(R_sc_km_vals)[:,2]
		self.vx_kms_vals = np.stack(R_sc_dot_kms_vals)[:,0]
		self.vy_kms_vals = np.stack(R_sc_dot_kms_vals)[:,1]
		self.vz_kms_vals = np.stack(R_sc_dot_kms_vals)[:,2]

		if test_file:
				self.run_test_file(test_file)
			

		if write_file:
				self.write_to_file(write_file)





	def et2spm(self, ET_vals):
			"""This converts ephemeris time to seconds past midnight.

			Args:
				ET_vals(arr):
					Array of ephemeris time.

			Returns:
				SPM_vals(arr):
					Array of seconds past midnight.
			"""
			SPM_vals = []
			for n in range(len(ET_vals)):
				ET = ET_vals[n]
				UTC_str = spice.et2utc(ET, 'ISOD', 16)
				hour = int(UTC_str[9:11])
				minute  = int(UTC_str[12:14])
				second = int(UTC_str[15:17])
				SPM = hour*3600. + minute*60. + second
				SPM_vals.append(SPM)
			return SPM_vals



	def spm2et(self, SPM_vals, DOY, year): 
		"""This converts seconds past midnight to ephemeris time/

		Args:
			SPM_vals(arr):
				Array of seconds past midnight.
			DOY(int):
				Day of year.
			year(int):
				Year.
		
		Returns:
			ET_vals(arr):
				Array of ephemeris times.
		"""

		ET_vals = []
		
		# Determine if leap year
		if year%4 == 0:
			days_per_year = 366
		else:
			days_per_year = 365

		for n in range(len(SPM_vals)):
			SPM = SPM_vals[n]

			this_hr = int(SPM/3600.)
			this_min = int((SPM - this_hr*3600.)/60.)
			this_sec = int(SPM - this_hr*3600. - this_min*60.)

			this_DOY = int(DOY + this_hr/24.)
			this_year = int(year + this_DOY/(days_per_year + 1))

			# Construct UTC string in ISOD format (ex: 1995-08T18:28:12)
			UTC_str = (str(this_year) + '-' + str(this_DOY) + 'T' +
					str(this_hr) + ':' + str(this_min) + ':' + str(this_sec))

			ET = spice.str2et(UTC_str)
			print(ET)
			ET_vals.append(ET)

		return ET_vals


	def get_pole(self, t_SET_et_vals):
		"""This calculates unit vector in pole direction from kernel constants.

		Args:
			t_SET_et_vals(arr):
				Array of spacecraft event time in ephemeris time.

		Returns:
			nhat_p_vals(arr):
				Array of unit vectors in pole direction/
		"""

		nhat_p_vals = []
		for n in t_SET_et_vals:
			bodynm = PLANET
			item = 'POLE_RA'
			maxn = 3
			dim1, pole_RA = spice.bodvrd(bodynm, item, maxn)

			bodynm = PLANET
			item = 'POLE_DEC'
			maxn = 3
			dim2, pole_DEC = spice.bodvrd(bodynm, item, maxn)

			dt_centuries = n / (spice.spd() * 365.25*100.)
			RAP = (pole_RA[0] + dt_centuries*pole_RA[1] + 
							dt_centuries**2*pole_RA[2])
			DEP = (pole_DEC[0] + dt_centuries*pole_DEC[1] + 
							dt_centuries**2*pole_DEC[2])

			inrange = 1.
			re = RAP * spice.rpd()
			dec = DEP * spice.rpd()
			nhat_p = spice.radrec(inrange, re, dec)
			nhat_p_vals.append(nhat_p)

		return nhat_p_vals

	def calc_dynam_frame(self, t_OET_spm_vals, year, DOY, DSN):
		"""This calculates ring intercept point using a dynamical frame.

		Args:
			t_OET_spm_vals(arr):
				Array of observed event times in seconds past midnight.
			year(int):
				Year of event.
			DOY(int):
				Day of year of event.
			DSN(str):
				ID of Deep Space Network receiver.

		Returns:
			rho_km_vals(arr):
				Array of ring intercept points in km.
			spoint_J2000_km_vals(arr):
			R_sc2dsn_km_vals(arr):
		"""
		
		rho_km_vals	= []
		spoint_J2000_km_vals = []
		R_sc2dsn_km_vals = []

		t_RET_et_vals = []
		t_SET_et_vals = []
	
		t_OET_et_vals = self.spm2et(t_OET_spm_vals, DOY, year)
		for n in range(len(t_OET_spm_vals)):
			t_OET_et = t_OET_et_vals[n]


			targ	= SPACECRAFT
			ref	= 'J2000'
			abcorr	= 'CN'
			obs	= DSN

			# spacecraft position relative to DSN station at
			# transmission time, and light travel time
			R_sc2dsn_km, ltime_sc2dsn = spice.spkpos(targ, t_OET_et, ref, 
							abcorr, obs)
			R_sc2dsn_km_vals.append(R_sc2dsn_km)
			
			t_SET_et = t_OET_et - ltime_sc2dsn
			t_SET_et_vals.append(t_SET_et)

			nhat_sc2dsn = spice.vhat(R_sc2dsn_km)

			# Retrieve radii from kernel pool
			body = SATURN_NAIF_ID
			item = 'RADII'
			dim = 3
			radii = spice.bodvcd(body, item, dim)
			
			# New radii to create ring plane frame
			new_radii = [1.e6, 1.e6, 1.e-5]

			# Replace old radii values with new
			name = SATURN_NAIF_RADII
			dvals = new_radii
			spice.pdpool(name, dvals)
			
			sincpt_method   = 'Ellipsoid'
			fixref          = 'IAU_SATURN'
			dref            = 'J2000'
			abcorr          = 'CN'
			obsrvr          = DSN
			dvec            = nhat_sc2dsn
			target 			= PLANET

            #  Compute intersection of signal with the ring plane frame
			spoint, t_RET_et, srfvec = spice.sincpt(sincpt_method, target,
					t_OET_et, fixref, abcorr, obsrvr, dref, dvec)
			
			t_RET_et_vals.append(t_RET_et)

            # Restore old valuoes of RADII to kernel pool
			spice.pdpool, 'RADII', radii

			# Convert ring plane intercept to J2000 frame
			frame_from = fixref
			frame_to = dref
			etfrom = t_RET_et
			etto = t_OET_et
			xform = spice.pxfrm2(frame_from, frame_to, etfrom, etto)
			
			spoint_J2000 = spice.mxv(xform, spoint)
			spoint_J2000_km_vals.append(spoint_J2000)

			rho_km = spice.vnorm(spoint_J2000)
			rho_km_vals.append(rho_km)



		self.OET_ET_vals = t_OET_et_vals
		self.RET_ET_vals = t_RET_et_vals
		self.SET_ET_vals = t_SET_et_vals


			
		return rho_km_vals, spoint_J2000_km_vals, R_sc2dsn_km_vals
	
	def calc_phi(self, spoint_J2000_km_vals, R_sc2dsn_km_vals):
		"""This calculates observed ring azimuth and ring longitude.

		Args:
			spoint_J2000_km_vals
			R_sc2dsn_km_vals

		Returns:
			phi_ORA_deg_vals
			phi_RL_deg_vals
		"""
		phi_RL_deg_vals = []
		phi_ORA_deg_vals = []
		nhat_p_vals = self.get_pole(self.SET_ET_vals)
		# transform intercept to ring plane frame
		for n in range(len(nhat_p_vals)):

			spoint_J2000_km = spoint_J2000_km_vals[n]
			R_dsn2sc_km = -1.*R_sc2dsn_km_vals[n]

			zaxis = [0., 0., 1.]
			axdef = nhat_p_vals[n]
			indexa = 3
			plndef = zaxis
			indexp = 2
			mm13 = spice.twovec(axdef, indexa, plndef, indexp)
			vec_RL = spice.mxv(mm13, spoint_J2000_km)
			
			radius_vec_RL, RA_vec_RL, DEC_vec_RL = spice.recrad(vec_RL)

			phi_RL_deg = RA_vec_RL * spice.dpr()
			phi_RL_deg_vals.append(phi_RL_deg)

			vec_ORA = spice.mxv(mm13, R_dsn2sc_km)
			vec_ORA[2] = 0.
			radius_vec_ORA, RA_vec_ORA, DEC_vec_ORA = spice.recrad(vec_ORA)
			phi_ORA_deg = (phi_RL_deg - RA_vec_ORA*spice.dpr() + 720.) % 360.
			phi_ORA_deg_vals.append(phi_ORA_deg)



		return phi_RL_deg_vals, phi_ORA_deg_vals

	def calc_D(self, t_RET_et_vals, t_SET_et_vals):
		"""This calculates distance from spacecraft to ring intercept point.

		Args:
			t_RET_et_vals
			t_SET_et_vals
		Returns:
			D_km_vals
		"""
		D_km_vals = []
		for n in range(len(t_RET_et_vals)):
				D_km = (t_RET_et_vals[n] - t_SET_et_vals[n]) * spice.clight()
				D_km_vals.append(D_km)

		return D_km_vals

	def calc_B(self, t_RET_et_vals, DSN):
		"""This calculates ring opening angle.

		Args:
			t_RET_et_vals
			DSN
		Returns:
			B_deg_vals
		"""
		ref = 'J2000'
		abcorr = 'CN'

		B_deg_vals = []
		nhat_p_vals = self.get_pole(self.SET_ET_vals)

		for n in range(len(t_RET_et_vals)):
			targ = DSN
			et = t_RET_et_vals[n]
			ref = 'J2000'
			abcorr = 'CN'
			obs = SPACECRAFT
			starg1, ltime1 = spice.spkpos(targ, et, ref, abcorr, obs)

			satpol = nhat_p_vals[n]

			
			v1 = starg1
			v2 = satpol
			B_rad = (np.pi/2.) - spice.vsep(v1, v2)
			B_deg_vals.append(B_rad * spice.dpr())

		return B_deg_vals

	def calc_F(self, D_km_vals, F_sky_hz_vals, B_deg_vals, phi_ORA_deg_vals):
		"""This calculates the Fresnel scale (see MTR1986 Eq 6).

		Args:
			D_km_vals
			F_sky_hz_vals
			B_deg_vals
			phi_ORA_deg_vals
		
		Returns:
			F_km_vals
		"""
		F_km_vals = []
		for n in range(len(D_km_vals)):
			lambda_sky = spice.clight() / F_sky_hz_vals[n]
			phi_ORA_rad = spice.rpd() * phi_ORA_deg_vals[n]
			B_rad = spice.rpd() * B_deg_vals[n]
			D_km = D_km_vals[n]

			F_km = np.sqrt((0.5* lambda_sky * D_km * 
							(1-(np.cos(B_rad))**2* 
							(np.sin(phi_ORA_rad))**2))/ 
							(np.sin(B_rad))**2)
			F_km_vals.append(F_km)
		return F_km_vals


	def calc_SC_state(self, t_SET_et_vals, DSN):
		"""This calculates spacecraft state vector in planetocentric frame.

		Args:
			t_SET_et_vals
			DSN

		Returns:
			R_sc_km_vals
			R_sc_dot_kms_vals
			R_imp_km_vals
		"""
		nhat_p_vals = self.get_pole(self.SET_ET_vals)
		R_sc_km_vals = []
		R_sc_dot_kms_vals = []
		R_imp_km_vals = []
		for n in range(len(nhat_p_vals)):
			t_SET_et = t_SET_et_vals[n]
			nhat_p = nhat_p_vals[n]
			# Saturn->Cassini state vector
			targ = SPACECRAFT
			et = t_SET_et
			ref = 'J2000'
			abcorr = 'NONE'
			obs = PLANET
			starg0, ltime0 = spice.spkezr(targ, et, ref, abcorr, obs)

			# Cassini->DSN state vector
			targ = DSN
			et = t_SET_et
			ref = 'J2000' 
			abcorr = 'XCN'
			obs = SPACECRAFT
			starg1, ltime1 = spice.spkezr(targ, et, ref, abcorr, obs)

			# Define planetocentric frame as z-axis in Saturn pole direction
			# and x-axis in Cassini to DSN direction
			zaxis = [0., 0., 1.]
			xaxis = [1., 0., 0.]

			axdef = nhat_p
			indexa = 3
			plndef = zaxis
			index = 2
			rotmat_z = spice.twovec(axdef, indexa, plndef, index)

			R_sat2sc_km = starg0[0:3]
			R_sc2dsn_km = starg1[0:3]
			R_sat2sc_dot_kms = starg0[3:6]
			

			R_sat2sc_km_z = spice.mxv(rotmat_z, R_sat2sc_km)
			R_sc2dsn_km_z = spice.mxv(rotmat_z, R_sc2dsn_km)
			R_sat2sc_dot_kms_z = spice.mxv(rotmat_z, R_sat2sc_dot_kms)

			nhat_sc2dsn = R_sc2dsn_km_z/np.linalg.norm(R_sc2dsn_km_z)

			nhat_sc2dsn[2] = 0.

			rot_angle = spice.vsep(nhat_sc2dsn, xaxis)

			rot_mat_sc2dsn = spice.rotate(rot_angle, 3)

			R_sat2sc_km_pcf = spice.mxv(rot_mat_sc2dsn, R_sat2sc_km_z)
			R_sc2dsn_km_pcf = spice.mxv(rot_mat_sc2dsn, R_sc2dsn_km_z)
			R_sat2sc_dot_kms_pcf = spice.mxv(rot_mat_sc2dsn, R_sat2sc_dot_kms_z)

			R_sc_km_vals.append(R_sat2sc_km_pcf)
			R_sc_dot_kms_vals.append(R_sat2sc_dot_kms_pcf)

			lindir = R_sc2dsn_km_pcf
			linpt = R_sat2sc_km_pcf
			point = [0., 0., 0.]

			pnear, distance = spice.nplnpt(linpt, lindir, point)
			R_imp_km_vals.append(distance)
		return R_sc_km_vals, R_sc_dot_kms_vals, R_imp_km_vals

	def run_test_file(self, test_file):
		"""This compares output results to a given test file.

		Args:
			test_file
		"""
#			EM_file = '/Volumes/jfong001/Research/TC2017/data/CORSS_8001_lien_resolution/EASYDATA/Rev07E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO.TAB'
		with open(test_file, 'r') as f:
			lines = f.readlines()
			t_OET_spm_vals_test = []
			t_RET_spm_vals_test = []
			t_SET_spm_vals_test = []
			rho_km_vals_test = []
			phi_RL_deg_vals_test = []
			phi_ORA_deg_vals_test = []
			B_deg_vals_test = []
			D_km_vals_test = []
			rho_dot_kms_vals_test = []
			phi_RL_dot_kms_vals_test = []
			F_km_vals_test = []
			R_imp_km_vals_test = []
			rx_km_vals_test = []
			ry_km_vals_test = []
			rz_km_vals_test = []
			vx_kms_vals_test = []
			vy_kms_vals_test = []
			vz_kms_vals_test = []
		
			t_OET_spm_vals_diff = []
			t_RET_spm_vals_diff = []
			t_SET_spm_vals_diff = []
			rho_km_vals_diff = []
			phi_RL_deg_vals_diff = []
			phi_ORA_deg_vals_diff = []
			B_deg_vals_diff = []
			D_km_vals_diff = []
			rho_dot_kms_vals_diff = []
			phi_RL_dot_kms_vals_diff = []
			F_km_vals_diff = []
			R_imp_km_vals_diff = []
			rx_km_vals_diff = []
			ry_km_vals_diff = []
			rz_km_vals_diff = []
			vx_kms_vals_diff = []
			vy_kms_vals_diff = []
			vz_kms_vals_diff = []

			delim = ','

			for line in lines:
				t_OET_spm_vals_test.append(float(line.split(delim)[0]))
				t_RET_spm_vals_test.append(float(line.split(delim)[1]))
				t_SET_spm_vals_test.append(float(line.split(delim)[2]))
				rho_km_vals_test.append(float(line.split(delim)[3]))
				phi_RL_deg_vals_test.append(float(line.split(delim)[4]))
				phi_ORA_deg_vals_test.append(float(line.split(delim)[5]))
				B_deg_vals_test.append(float(line.split(delim)[6]))
				D_km_vals_test.append(float(line.split(delim)[7]))
				rho_dot_kms_vals_test.append(float(line.split(delim)[8]))
				phi_RL_dot_kms_vals_test.append(float(line.split(delim)[9]))
				F_km_vals_test.append(float(line.split(delim)[10]))
				R_imp_km_vals_test.append(float(line.split(delim)[11]))
				rx_km_vals_test.append(float(line.split(delim)[12]))
				ry_km_vals_test.append(float(line.split(delim)[13]))
				rz_km_vals_test.append(float(line.split(delim)[14]))
				vx_kms_vals_test.append(float(line.split(delim)[15]))
				vy_kms_vals_test.append(float(line.split(delim)[16]))
				vz_kms_vals_test.append(float(line.split(delim)[17]))

			for n in range(len(self.t_OET_spm_vals)):
				t_OET_spm_vals_diff.append(self.t_OET_spm_vals[n]-
								t_OET_spm_vals_test[n])
				t_RET_spm_vals_diff.append(self.t_RET_SPM_vals[n] - 
								t_RET_spm_vals_test[n])
				t_SET_spm_vals_diff.append(self.t_SET_SPM_vals[n] - 
								t_SET_spm_vals_test[n])
				rho_km_vals_diff.append(self.rho_km_vals[n] - 
								rho_km_vals_test[n])
				phi_RL_deg_vals_diff.append(self.phi_RL_deg_vals[n] - 
								phi_RL_deg_vals_test[n])
				phi_ORA_deg_vals_diff.append(self.phi_ORA_deg_vals[n] - 
								phi_ORA_deg_vals_test[n])
				B_deg_vals_diff.append(self.B_deg_vals[n] - 
								B_deg_vals_test[n])
				D_km_vals_diff.append(self.D_km_vals[n] - 
								D_km_vals_test[n])
				rho_dot_kms_vals_diff.append(self.rho_dot_kms_vals[n] - 
								rho_dot_kms_vals_test[n])
				phi_RL_dot_kms_vals_diff.append(self.phi_RL_dot_kms_vals[n] - 
								phi_RL_dot_kms_vals_test[n])
				F_km_vals_diff.append(self.F_km_vals[n] - 
								F_km_vals_test[n])
				R_imp_km_vals_diff.append(self.R_imp_km_vals[n] - 
								R_imp_km_vals_test[n])
				rx_km_vals_diff.append(self.rx_km_vals[n] - 
								rx_km_vals_test[n])
				ry_km_vals_diff.append(self.ry_km_vals[n] - 
								ry_km_vals_test[n])
				rz_km_vals_diff.append(self.rz_km_vals[n] - 
								rz_km_vals_test[n])
				vx_kms_vals_diff.append(self.vx_kms_vals[n] - 
								vx_kms_vals_test[n])
				vy_kms_vals_diff.append(self.vy_kms_vals[n] - 
								vy_kms_vals_test[n])
				vz_kms_vals_diff.append(self.vz_kms_vals[n] - 
								vz_kms_vals_test[n])
			print('MINIMUM AND MAXIMUM DIFFERENCES:')
			print('t_OET_spm: ', min(t_OET_spm_vals_diff),
							max(t_OET_spm_vals_diff))
			print('t_RET_spm: ', min(t_RET_spm_vals_diff), 
							max(t_RET_spm_vals_diff))
			print('t_SET_spm: ', min(t_SET_spm_vals_diff), 
							max(t_SET_spm_vals_diff))
			print('rho_km: ', min(rho_km_vals_diff), max(rho_km_vals_diff))
			print('phi_RL_deg: ', min(phi_RL_deg_vals_diff), 
							max(phi_RL_deg_vals_diff))
			print('phi_ORA_deg: ', min(phi_ORA_deg_vals_diff), 
							max(phi_ORA_deg_vals_diff))
			print('B_deg: ', min(B_deg_vals_diff), max(B_deg_vals_diff))
			print('D_km: ', min(D_km_vals_diff), max(D_km_vals_diff))
			print('rho_dot_kms: ', min(rho_dot_kms_vals_diff), 
							max(rho_dot_kms_vals_diff))
			print('phi_RL_dot_kms: ', min(phi_RL_dot_kms_vals_diff), 
							max(phi_RL_dot_kms_vals_diff))
			print('F_km: ', min(F_km_vals_diff), max(F_km_vals_diff))
			print('R_imp_km: ', min(R_imp_km_vals_diff), \
							max(R_imp_km_vals_diff))
			print('rx_km: ', min(rx_km_vals_diff), max(rx_km_vals_diff))
			print('ry_km: ', min(ry_km_vals_diff), max(ry_km_vals_diff))
			print('rz_km: ', min(rz_km_vals_diff), max(rz_km_vals_diff))
			print('vx_kms: ', min(vx_kms_vals_diff), max(vx_kms_vals_diff))
			print('vy_kms: ', min(vy_kms_vals_diff), max(vy_kms_vals_diff))
			print('vz_kms: ', min(vz_kms_vals_diff), max(vz_kms_vals_diff))

			# Plot fractional error
			t_OET_spm_vals_err = [a/b for a,b in zip(t_OET_spm_vals_diff,
					t_OET_spm_vals_test)]
			t_RET_spm_vals_err = [a/b for a,b in zip(t_RET_spm_vals_diff,
					t_RET_spm_vals_test)]
			t_SET_spm_vals_err = [a/b for a,b in zip(t_SET_spm_vals_diff,
					t_SET_spm_vals_test)]
			rho_km_vals_err = [a/b for a,b in zip(rho_km_vals_diff, 
					rho_km_vals_test)]
			phi_RL_deg_vals_err = [a/b for a,b in zip(phi_RL_deg_vals_diff,
					phi_RL_deg_vals_test)]
			phi_ORA_deg_vals_err = [a/b for a,b in zip(phi_ORA_deg_vals_diff,
					phi_ORA_deg_vals_test)]
			B_deg_vals_err = [a/b for a,b in zip(B_deg_vals_diff,
					B_deg_vals_test)]
			D_km_vals_err = [a/b for a,b in zip(D_km_vals_diff,
					D_km_vals_test)]
			rho_dot_kms_vals_err = [a/b for a,b in zip(rho_dot_kms_vals_diff,
					rho_dot_kms_vals_test)]
			phi_RL_dot_kms_vals_err = [a/b for a,b in zip(
					phi_RL_dot_kms_vals_diff, phi_RL_dot_kms_vals_test)]
			F_km_vals_err = [a/b for a,b in zip(F_km_vals_diff,
					F_km_vals_test)]
			R_imp_km_vals_err = [a/b for a,b in zip(R_imp_km_vals_diff,
					R_imp_km_vals_test)]
			rx_km_vals_err = [a/b for a,b in zip(rx_km_vals_diff, 
					rx_km_vals_test)]
			ry_km_vals_err = [a/b for a,b in zip(ry_km_vals_diff,
					ry_km_vals_test)]
			rz_km_vals_err = [a/b for a,b in zip(rz_km_vals_diff,
					rz_km_vals_test)]
			vx_kms_vals_err = [a/b for a,b in zip(vx_kms_vals_diff,
					vx_kms_vals_test)]
			vy_kms_vals_err = [a/b for a,b in zip(vy_kms_vals_diff,
					vy_kms_vals_test)]
			vz_kms_vals_err = [a/b for a,b in zip(vz_kms_vals_diff,
					vz_kms_vals_test)]
			fig1, axes1 = plt.subplots(nrows=10, ncols=1, sharex=True)
			x = self.t_OET_spm_vals # for easier plotting...
			axes1[0].set_title('Fractional Error')
			axes1[0].plot(x, t_OET_spm_vals_err)
			axes1[1].plot(x, t_RET_spm_vals_err)
			axes1[2].plot(x, t_SET_spm_vals_err)
			axes1[3].plot(x, rho_km_vals_err)
			axes1[4].plot(x, phi_RL_deg_vals_err)
			axes1[5].plot(x, phi_ORA_deg_vals_err)
			axes1[6].plot(x, B_deg_vals_err)
			axes1[7].plot(x, D_km_vals_err)
			axes1[8].plot(x, rho_dot_kms_vals_err)
			axes1[9].plot(x, phi_RL_dot_kms_vals_err)
			axes1[9].set_xlabel('Observed Event Time (1000 s)')

			fig2, axes2 = plt.subplots(nrows=8, ncols=1, sharex=True)
			axes2[0].set_title('Fractional Error')
			axes2[0].plot(x, F_km_vals_err)
			axes2[1].plot(x, R_imp_km_vals_err)
			axes2[2].plot(x, rx_km_vals_err)
			axes2[3].plot(x, ry_km_vals_err)
			axes2[4].plot(x, rz_km_vals_err)
			axes2[5].plot(x, vx_kms_vals_err)
			axes2[6].plot(x, vy_kms_vals_err)
			axes2[7].plot(x, vz_kms_vals_err)
			plt.show()


	def write_to_file(self, write_file):
		"""This writes output results to a text file.

		Args:
			write_file
		"""
		format_str = ['%20.9E'] * 18
		input_text = np.c_[self.t_OET_spm_vals, 
						 self.t_RET_SPM_vals, 
						 self.t_SET_SPM_vals, 
						 self.rho_km_vals, 
						 self.phi_RL_deg_vals, 
						 self.phi_ORA_deg_vals,
						 self.B_deg_vals, 
						 self.D_km_vals, 
						 self.rho_dot_kms_vals, 
						 self.phi_RL_dot_kms_vals, 
						 self.F_km_vals, 
						 self.R_imp_km_vals, 
						 self.rx_km_vals, 
						 self.ry_km_vals, 
						 self.rz_km_vals, 
						 self.vx_kms_vals, 
						 self.vy_kms_vals, 
						 self.vz_kms_vals]
		np.savetxt(write_file, input_text, delimiter = ',', fmt=format_str) 


def main():
	start_time = time.time()
	RSRfile = '/Volumes/jfong001/Research/TC2017/data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'

	kernel_dir = '/Volumes/jfong001/Research/TC2017/kernels/'
	kernels = ['050606R_SCPSE_05114_05132.bsp'         
			,'cpck26Feb2009.tpc'                    
			,'naif0009.tls'				
			,'earthstns_itrf93_050714.bsp'          
			,'de421.bsp'                            
			,'earth_000101_090604_090313.bpc']
	kernel_list = map(lambda x: kernel_dir + x, kernels)

	rsr_reader = RSRReader(RSRfile)
	rsr_reader.read_hdr()
#	print(rsr_reader.__dict__)
	rsr_reader.get_f_sky_pred()

#	test_file = '/Volumes/jfong001/Research/TC2017/data/CORSS_8001_lien_resolution/EASYDATA/Rev07E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO.TAB'
	write_file = 'rev7E_geometry_test_file.tab'
	test_file = '/Volumes/jfong001/Research/TC2017/jfong/progress_report/Rev007_E_X43_GEO_v2_E2012.tab'
	write_file = 0
	test_file = 0


	geom = Geom(rsr_reader, kernel_list, test_file=test_file, write_file=
			write_file)
	end_time = time.time()
	print('Run time: ', end_time-start_time)
if __name__ == '__main__':
	main()


