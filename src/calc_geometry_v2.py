
'''

calc_geometry_v2.py
Purpose: Class Geom whose attributes are the occultation geometry parameters.
Revisions:
	2018 Feb 12 - jfong - original
	2018 Feb 17 - jfong - compatible with jwf_read_rsr_v2
	2018 Feb 22 - jfong - add test_file kwd that compares results to given file
	2018 Mar 02 - jfong - test_file kwd plots fractional error
	2018 Mar 13 - jfong - v2 - use gjs_spm_to_et.py
							 - use numpy arrays instead of lists
							 - remove test_file and write_file kwds
	2018 Mar 29 - jfong - use et_to_spm.py
						- use lists for 2d numpy arrays to avoid strided arrays
							- spoint_J2000_km, R_sc2dsn_km, nhat_p
							  R_sc_km_vals, R_sc_dot_kms_vals, R_imp_km_vals

'''

import spiceypy as spice
import numpy as np
from rsr_reader import RSRReader
from spm_to_et import spm_to_et
from et_to_spm import et_to_spm
from make_geo_file import make_geo_file
#from jwf_compare_geo_file import compare_geo_file
#import matplotlib.pyplot as plt
import time

SPACECRAFT              = 'Cassini'
PLANET                  = 'Saturn'
SATURN_NAIF_ID          = 699
SATURN_NAIF_RADII       = 'BODY699_RADII'

class Geometry(object):

	"""This is an object that calculates occultation geometry from event time.
	"""


	def __init__(self, rsr_header, kernels, pt_per_sec=None):
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
		self.kernels = kernels

		self.t_oet_spm_vals			= None
		self.t_ret_et_vals			= None
		self.t_set_et_vals			= None
		self.rho_km_vals			= None
		self.phi_rl_deg_vals		= None
		self.phi_ora_deg_vals		= None
		self.B_deg_vals				= None
		self.D_km_vals				= None
		self.rho_dot_kms_vals		= None
		self.phi_rl_dot_kms_vals	= None
		self.F_km_vals				= None
		self.R_imp_km_vals			= None
		self.rx_km_vals				= None
		self.ry_km_vals				= None
		self.rz_km_vals				= None
		self.vx_kms_vals			= None
		self.vy_kms_vals			= None
		self.vz_kms_vals			= None

		self.__pole					= None

		if pt_per_sec == None:
			pt_per_sec = 1.

		year	= rsr_header.year
		doy		= rsr_header.doy
		dsn		= rsr_header.dsn
		band	= rsr_header.band
		
		# Create new spm array with defined points per second		
		spm_start = rsr_header.spm_vals[0]
		spm_end = rsr_header.spm_vals[-1]
		step = 1./pt_per_sec
		t_oet_spm_vals = np.arange(spm_start, spm_end, step)

		(f_spm, f_sky) = rsr_header.get_f_sky_pred()
		F_sky_hz_vals = np.interp(t_oet_spm_vals, f_spm, f_sky)


		# Load kernels
		self.__load_kernels(kernels)


		rho_km_vals, spoint_J2000_km_vals, R_sc2dsn_km_vals = \
				self.calc_dynam_frame(t_oet_spm_vals, year, 
				doy, dsn, kernels=kernels)

		self.__pole = self.__get_pole(self.t_set_et_vals)

		phi_rl_deg_vals, phi_ora_deg_vals = self.calc_phi(self.t_set_et_vals,
				spoint_J2000_km_vals, R_sc2dsn_km_vals)

		D_km_vals = self.calc_D(self.t_ret_et_vals, self.t_set_et_vals)

		B_deg_vals = self.calc_B(self.t_set_et_vals, dsn)

		F_km_vals = self.calc_F(D_km_vals, F_sky_hz_vals, B_deg_vals, 
				phi_ora_deg_vals)

		R_sc_km_vals, R_sc_dot_kms_vals, R_imp_km_vals = \
				self.calc_sc_state(self.t_set_et_vals, dsn)

		t_ret_spm_vals = et_to_spm(self.t_ret_et_vals)
		t_set_spm_vals = et_to_spm(self.t_set_et_vals)


		rho_dot_kms_vals = rho_km_vals - np.roll(rho_km_vals, 1)
		phi_rl_rad_vals = [x*spice.rpd() for x in phi_rl_deg_vals]
		phi_rl_dot_rads_vals = phi_rl_rad_vals - np.roll(phi_rl_rad_vals, 1)
		phi_rl_dot_kms_vals	= [a*b for a,b in zip(rho_km_vals, 
					phi_rl_dot_rads_vals)]

		
		self.t_oet_spm_vals = t_oet_spm_vals
		self.t_ret_spm_vals = t_ret_spm_vals
		self.t_set_spm_vals = t_set_spm_vals
		self.rho_km_vals = rho_km_vals
		self.phi_rl_deg_vals = phi_rl_deg_vals
		self.phi_ora_deg_vals = phi_ora_deg_vals
		self.D_km_vals = D_km_vals
		self.B_deg_vals = B_deg_vals
		self.rho_dot_kms_vals = rho_dot_kms_vals
		self.phi_rl_dot_kms_vals = phi_rl_dot_kms_vals
		self.F_km_vals = F_km_vals
		self.R_imp_km_vals = R_imp_km_vals
		self.rx_km_vals = np.stack(R_sc_km_vals)[:,0]
		self.ry_km_vals = np.stack(R_sc_km_vals)[:,1]
		self.rz_km_vals = np.stack(R_sc_km_vals)[:,2]
		self.vx_kms_vals = np.stack(R_sc_dot_kms_vals)[:,0]
		self.vy_kms_vals = np.stack(R_sc_dot_kms_vals)[:,1]
		self.vz_kms_vals = np.stack(R_sc_dot_kms_vals)[:,2]

	def __load_kernels(self, kernels):
		spice.kclear()
		for kernel in kernels:
			spice.furnsh(kernel)

	def __get_pole(self, t_set_et_vals, kernels=None):
		"""This calculates unit vector in pole direction from kernel constants.

		Args:
			t_set_et_vals(arr):
				Array of spacecraft event time in ephemeris time.

		Returns:
			nhat_p_vals(arr):
				Array of unit vectors in pole direction/
		"""
		if kernels:
			self.__load_kernels(kernels)

		npts = len(t_set_et_vals)
		nhat_p_vals = []
		for n in range(npts):

			t_set_et = t_set_et_vals[n]

			bodynm = PLANET
			item = 'POLE_RA'
			maxn = 3
			dim1, pole_RA = spice.bodvrd(bodynm, item, maxn)

			bodynm = PLANET
			item = 'POLE_DEC'
			maxn = 3
			dim2, pole_DEC = spice.bodvrd(bodynm, item, maxn)

			dt_centuries = t_set_et / (spice.spd() * 365.25*100.)
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

	def calc_dynam_frame(self, t_oet_spm_vals, year, doy, dsn, kernels=None):
		"""This calculates ring intercept point using a dynamical frame.

		Args:
			t_oet_spm_vals(arr):
				Array of observed event times in seconds past midnight.
			year(int):
				Year of event.
			doy(int):
				Day of year of event.
			dsn(str):
				ID of Deep Space Network receiver.

		Returns:
			rho_km_vals(arr):
				Array of ring intercept points in km.
			spoint_J2000_km_vals(arr):
			R_sc2dsn_km_vals(arr):
		"""
		npts = len(t_oet_spm_vals)
		rho_km_vals	= np.zeros(npts)
		spoint_J2000_km_vals = [] #np.zeros((3, npts))
		R_sc2dsn_km_vals = [] # np.zeros((3, npts))

		t_ret_et_vals = np.zeros(npts)
		t_set_et_vals = np.zeros(npts)
	
		t_oet_et_vals = spm_to_et(t_oet_spm_vals, doy, year, kernels=kernels)
		for n in range(npts):
			t_oet_et = t_oet_et_vals[n]


			targ	= SPACECRAFT
			ref	= 'J2000'
			abcorr	= 'CN'
			obs	= dsn

			# spacecraft position relative to dsn station at
			# transmission time, and light travel time
			R_sc2dsn_km, ltime_sc2dsn = spice.spkpos(targ, t_oet_et, ref, 
							abcorr, obs)
		
			R_sc2dsn_km_vals.append(R_sc2dsn_km)
			
			t_set_et = t_oet_et - ltime_sc2dsn
			t_set_et_vals[n] = t_set_et

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
			
            #  Compute intersection of signal with the ring plane frame
			sincpt_method   = 'Ellipsoid'
			fixref          = 'IAU_SATURN'
			dref            = 'J2000'
			abcorr          = 'CN'
			obsrvr          = dsn
			dvec            = nhat_sc2dsn
			target 			= PLANET
			spoint, t_ret_et, srfvec = spice.sincpt(sincpt_method, target,
					t_oet_et, fixref, abcorr, obsrvr, dref, dvec)
			
			t_ret_et_vals[n] = t_ret_et

            # Restore old valuoes of RADII to kernel pool
			spice.pdpool, 'RADII', radii

			# Convert ring plane intercept to J2000 frame
			frame_from = fixref
			frame_to = dref
			etfrom = t_ret_et
			etto = t_oet_et
			xform = spice.pxfrm2(frame_from, frame_to, etfrom, etto)
			
			spoint_J2000_km = spice.mxv(xform, spoint)
#			spoint_J2000_km_vals[0:3, n] = spoint_J2000_km
			spoint_J2000_km_vals.append(spoint_J2000_km)

			rho_km = spice.vnorm(spoint_J2000_km)
			rho_km_vals[n] = rho_km



		self.t_oet_et_vals = t_oet_et_vals
		self.t_ret_et_vals = t_ret_et_vals
		self.t_set_et_vals = t_set_et_vals


			
		return rho_km_vals, spoint_J2000_km_vals, R_sc2dsn_km_vals
	
	def calc_phi(self, t_set_et_vals, spoint_J2000_km_vals, R_sc2dsn_km_vals):
		"""This calculates observed ring azimuth and ring longitude.

		Args:
			spoint_J2000_km_vals
			R_sc2dsn_km_vals

		Returns:
			phi_ora_deg_vals
			phi_rl_deg_vals
		"""
		npts = len(t_set_et_vals)
		phi_rl_deg_vals = np.zeros(npts)
		phi_ora_deg_vals = np.zeros(npts)

		nhat_p_vals = self.__pole
		# transform intercept to ring plane frame
		for n in range(npts):

			spoint_J2000_km = spoint_J2000_km_vals[n]
			R_dsn2sc_km = -1.*R_sc2dsn_km_vals[n]

#			nhat_p = np.reshape(nhat_p_vals[0:3, n], (1,3))
			nhat_p = nhat_p_vals[n]

			zaxis = [0., 0., 1.]
			axdef = nhat_p
			indexa = 3
			plndef = zaxis
			indexp = 2
			mm13 = spice.twovec(axdef, indexa, plndef, indexp)
			vec_RL = spice.mxv(mm13, spoint_J2000_km)
			
			radius_vec_RL, RA_vec_RL, DEC_vec_RL = spice.recrad(vec_RL)

			phi_rl_deg = RA_vec_RL * spice.dpr()
			phi_rl_deg_vals[n] = phi_rl_deg

			vec_ORA = spice.mxv(mm13, R_dsn2sc_km)
			vec_ORA[2] = 0.
			radius_vec_ORA, RA_vec_ORA, DEC_vec_ORA = spice.recrad(vec_ORA)
			phi_ora_deg = (phi_rl_deg - RA_vec_ORA*spice.dpr() + 720.) % 360.
			phi_ora_deg_vals[n] = phi_ora_deg




		return phi_rl_deg_vals, phi_ora_deg_vals

	def calc_D(self, t_ret_et_vals, t_set_et_vals):
		"""This calculates distance from spacecraft to ring intercept point.

		Args:
			t_ret_et_vals
			t_set_et_vals
		Returns:
			D_km_vals
		"""
		npts = len(t_ret_et_vals)
		D_km_vals = np.zeros(npts)
		for n in range(npts):
			D_km = (t_ret_et_vals[n] - t_set_et_vals[n]) * spice.clight()
			D_km_vals[n] = D_km

		return D_km_vals

	def calc_B(self, t_ret_et_vals, dsn):
		"""This calculates ring opening angle.

		Args:
			t_ret_et_vals
			dsn
		Returns:
			B_deg_vals
		"""
		ref = 'J2000'
		abcorr = 'CN'
		
		npts = len(t_ret_et_vals)

		B_deg_vals = np.zeros(npts)
		nhat_p_vals = self.__pole

		for n in range(npts):
			targ = dsn
			et = t_ret_et_vals[n]
			ref = 'J2000'
			abcorr = 'CN'
			obs = SPACECRAFT
			starg1, ltime1 = spice.spkpos(targ, et, ref, abcorr, obs)

			satpol = nhat_p_vals[n]

			
			v1 = starg1
			v2 = satpol
			B_rad = (np.pi/2.) - spice.vsep(v1, v2)
			B_deg_vals[n] = B_rad * spice.dpr()

		return B_deg_vals

	def calc_F(self, D_km_vals, F_sky_hz_vals, B_deg_vals, phi_ora_deg_vals):
		"""This calculates the Fresnel scale (see MTR1986 Eq 6).

		Args:
			D_km_vals
			F_sky_hz_vals
			B_deg_vals
			phi_ora_deg_vals
		
		Returns:
			F_km_vals
		"""
		npts = len(D_km_vals)
		F_km_vals = np.zeros(npts)

		for n in range(npts):
			lambda_sky = spice.clight() / F_sky_hz_vals[n]
			phi_ORA_rad = spice.rpd() * phi_ora_deg_vals[n]
			B_rad = spice.rpd() * B_deg_vals[n]
			D_km = D_km_vals[n]

			F_km = np.sqrt((0.5* lambda_sky * D_km * 
							(1-(np.cos(B_rad))**2* 
							(np.sin(phi_ORA_rad))**2))/ 
							(np.sin(B_rad))**2)
			F_km_vals[n] = F_km
		return F_km_vals


	def calc_sc_state(self, t_set_et_vals, dsn):
		"""This calculates spacecraft state vector in planetocentric frame.

		Args:
			t_set_et_vals
			dsn

		Returns:
			R_sc_km_vals
			R_sc_dot_kms_vals
			R_imp_km_vals
		"""
		npts = len(t_set_et_vals)
		nhat_p_vals = self.__pole
		R_sc_km_vals = []
		R_sc_dot_kms_vals = []
		R_imp_km_vals = []
		for n in range(npts):
			t_set_et = t_set_et_vals[n]
			nhat_p = nhat_p_vals[n]
			# Saturn->Cassini state vector
			targ = SPACECRAFT
			et = t_set_et
			ref = 'J2000'
			abcorr = 'NONE'
			obs = PLANET
			starg0, ltime0 = spice.spkezr(targ, et, ref, abcorr, obs)

			# Cassini->dsn state vector
			targ = dsn
			et = t_set_et
			ref = 'J2000' 
			abcorr = 'XCN'
			obs = SPACECRAFT
			starg1, ltime1 = spice.spkezr(targ, et, ref, abcorr, obs)

			# Define planetocentric frame as z-axis in Saturn pole direction
			# and x-axis in Cassini to dsn direction
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


def main():
	start_time = time.time()
	RSRfile = '/Volumes/jfong001/Research/TC2017/data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'

	geo_file = 'jwf_calc_geometry_v2_output.tab'

	kernels_dir = '/Volumes/jfong001/Research/TC2017/kernels~/'
	kernels_files = ['050606R_SCPSE_05114_05132.bsp'         
			,'cpck26Feb2009.tpc'                    
			,'naif0009.tls'				
			,'earthstns_itrf93_050714.bsp'          
			,'de421.bsp'                            
			,'earth_000101_090604_090313.bpc']
#	kernel_list = map(lambda x: kernel_dir + x, kernels)
	kernels = [kernels_dir + this_kernel for this_kernel in kernels_files]

	rsr_reader = RSRReader(RSRfile)

	test_file = '/Volumes/jfong001/Research/TC2017/jfong/progress_report/Rev007_E_X43_GEO_v2_E2012.tab'
	test_file = 0


	geom = Geometry(rsr_reader, kernels)
	#print(geom.__dict__)
	make_geo_file(geom, geo_file)
	end_time = time.time()
	print('Run time: ', end_time-start_time)
if __name__ == '__main__':
	main()


