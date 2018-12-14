"""

Purpose:
         Create a class whose attributes have all the necessary inputs for
         performing a Fresnel inversion usng DiffractionCorrection, given
         instances of the classes RSRReader, Geometry, and Calibration.

"""

import copy
import numpy as np
import pdb
from scipy.interpolate import splrep, splev, interp1d
import sys

from ..rsr_reader.rsr_reader import RSRReader
from .resample_IQ import resample_IQ
from .calc_tau_thresh import calc_tau_thresh
from ..tools.history import write_history_dict
from ..tools.write_output_files import write_output_files

sys.path.append('../..')
import rss_ringoccs as rss
sys.path.remove('../..')

class DiffractionLimitedProfile(object):
    """
    Framework for an object class containing the diffraction-limited optical depth
    profile (DLP) and related attributes.

    Arguments:
        :rsr_inst (*object*): Instance of RSRReader class
        :geo_inst (*object*): Instance of Geometry class
        :cal_inst (*object*): Instance of Calibration class
        :dr_km (*float*): radial sampling rate :math:`\\Delta\\rho` for the DLP
                            in kilometers. DLP radial *resolution* is the Nyquist
                            radial sampling, i.e., twice the input value of `dr_km`,
                            meaning that this will affect the minimum resolution
                            of the diffraction-reconstructed profile. Value for
                            `dr_km` can range from 0.05 km to 0.75 km for the
                            reconstruction resolutions supported by `rss_ringoccs`.
                            PDS sampling rate is 0.25 km, which gives a DLP resolution of 0.5 km.

    Keyword Arguments:
        :dr_km_tol (*float*): Maximum distance from an integer number of
                    dr_desired that the first rho value will be at. Makes the final
                    set of rho values more regular-looking. For example, if you say
                    dr_km_tol=0.01 with a dr_desired of 0.25, your final set of rho
                    values might look something like [70000.26, 70000.51, ...]
        :verbose (*bool*): When True, turns on verbose output. Default is False.
        :write_file (*bool*): When True, writes processing results to file. Default is True.
        :profile_range (*list*): 1x2 list specifying the radial limits in km of on the
                    occultation. Default is [65000,150000].

    Attributes:
        :snr (*np.ndarray*): Signal-to-noise ratio
        :tau_thresh (*np.ndarray*): threshold optical depth
        :spm_thresh (*np.ndarray*): SPM at which threshold optical depth is computed
        :rho_thresh (*np.ndarray*): radius at which threshold optical depth is computed
        :dr_km (*float*): raw DLP sampling rate
    """
    def __init__(self, rsr_inst, geo_inst, cal_inst, dr_km,
            dr_km_tol=0.01, verbose=False, write_file=True,
            profile_range=[65000., 150000.]):


        # Extract necessary information from instance inputs
        spm_full = rsr_inst.spm_vals
        IQ_m = rsr_inst.IQ_m

        spm_geo = geo_inst.t_oet_spm_vals
        rho_km_geo = geo_inst.rho_km_vals
        rho_dot_kms_geo = geo_inst.rho_dot_kms_vals
        prof_dir = geo_inst.rev_info['prof_dir']
        B_deg_geo = geo_inst.B_deg_vals
        F_km_geo = geo_inst.F_km_vals
        t_ret_spm_geo = geo_inst.t_ret_spm_vals
        t_set_spm_geo = geo_inst.t_set_spm_vals
        phi_ora_deg_geo = geo_inst.phi_ora_deg_vals
        phi_rl_deg_geo = geo_inst.phi_rl_deg_vals
        D_km_geo = geo_inst.D_km_vals

        spm_cal = cal_inst.t_oet_spm_vals
        f_sky_pred_cal = cal_inst.f_sky_hz_vals
        p_free_cal = cal_inst.p_free_vals
        IQ_c = cal_inst.IQ_c




        # interpolate coef to convert SPM to rho
        spm_to_rho = splrep(spm_geo, rho_km_geo)
        rho_full = splev(spm_full, spm_to_rho)

        # get predicted sky frequency for cal SPM
        dummy_spm, f_sky_pred_file = rsr_inst.get_f_sky_pred(f_spm=spm_cal)
        f_offset_fit_cal = f_sky_pred_cal - f_sky_pred_file
        rho_km_cal = splev(spm_cal, spm_to_rho)

        # resample IQ_c, but not to exact rho boundaries
        rho_km_desired, IQ_c_desired = resample_IQ(rho_full, IQ_c,
                dr_km, dr_km_tol=dr_km_tol)

        # Slice rho_km_desired to profile_range values
        ind1 = list(rho_km_desired).index(min(list(rho_km_desired),
                    key=lambda x:abs(x-profile_range[0])))

        ind2 = list(rho_km_desired).index(min(list(rho_km_desired),
                    key=lambda x:abs(x-profile_range[1])))

        rho_km_desired = rho_km_desired[ind1:ind2+1]
        IQ_c_desired = IQ_c_desired[ind1:ind2+1]

        # interpolate coef to convert rho to SPm
        spm_desired_func = interp1d(rho_km_geo, spm_geo,
                fill_value='extrapolate')
        spm_desired = spm_desired_func(rho_km_desired)

        p_free_interp_func = interp1d(rho_km_cal, p_free_cal,
                fill_value='extrapolate')
        p_free = p_free_interp_func(rho_km_desired)

        p_norm_vals = (abs(IQ_c_desired)**2)/p_free

        phase_rad_vals = np.arctan2(np.imag(IQ_c_desired),
                np.real(IQ_c_desired))

        # compute threshold optical depth
        tau_thresh_inst = calc_tau_thresh(rsr_inst,geo_inst,cal_inst,res_km=0.5)
        self.snr = tau_thresh_inst.snr
        self.tau_thresh = tau_thresh_inst.tau_thresh
        self.spm_thresh = tau_thresh_inst.spm_vals
        self.rho_thresh = tau_thresh_inst.rho_vals

        self.interp_and_set_attr(rho_km_desired, spm_desired, p_norm_vals,
                spm_cal, phase_rad_vals, spm_geo, rho_dot_kms_geo,
                B_deg_geo, F_km_geo, t_ret_spm_geo, t_set_spm_geo,
                D_km_geo, phi_ora_deg_geo, phi_rl_deg_geo, f_sky_pred_cal)


        # if set, write output data and label file
        self.rev_info = geo_inst.rev_info
        self.dr_km = dr_km


        input_vars = {
                'rsr_inst': rsr_inst.history,
                'geo_inst': geo_inst.history,
                'cal_inst': cal_inst.history,
                'dr_km': dr_km}

        input_kwds = {
                'dr_km_tol': dr_km_tol,
                'profile_range': profile_range}

        self.history = write_history_dict(input_vars, input_kwds, __file__)
        if write_file:
            write_output_files(self)



    def interp_and_set_attr(self, rho_km_desired, spm_desired,
            p_norm_vals, spm_cal, phase_rad_vals, spm_geo, rho_dot_kms_geo,
            B_deg_vals, F_km_geo, t_ret_geo, t_set_geo,
            D_km_geo, phi_ora_deg_vals, phi_rl_deg_vals, f_sky_pred_cal):

        B_rad_geo = np.radians(B_deg_vals)
        phi_ora_rad_geo = np.radians(phi_ora_deg_vals)
        phi_rl_rad_geo = np.radians(phi_rl_deg_vals)

        # ring opening angle at final spacing
        spm_to_B = splrep(spm_geo, B_rad_geo)
        B_rad_vals_interp = splev(spm_desired, spm_to_B)

        # spacecraft - rip distance at final spacing
        spm_to_D = splrep(spm_geo, D_km_geo)
        D_km_vals_interp = splev(spm_desired, spm_to_D)

        # Fresnel scale at final spacing
        spm_to_F = splrep(spm_geo, F_km_geo)
        F_km_vals_interp = splev(spm_desired, spm_to_F)

        # sky frequency at final spacing
        spm_to_fsky = splrep(spm_cal, f_sky_pred_cal)
        f_sky_hz_vals_interp = splev(spm_desired, spm_to_fsky)

        # obvserved ring azimuth at final spacing
        spm_to_phi_ora = splrep(spm_geo, phi_ora_rad_geo)
        phi_ora_rad_vals_interp = splev(spm_desired, spm_to_phi_ora)

        # obvserved ring azimuth at final spacing
        spm_to_phi_rl = splrep(spm_geo, phi_rl_rad_geo)
        phi_rl_rad_vals_interp = splev(spm_desired, spm_to_phi_rl)

        # ring event time at final spacing
        spm_to_ret = splrep(spm_geo, t_ret_geo)
        t_ret_spm_vals_interp = splev(spm_desired, spm_to_ret)

        # spacecraft event time at final spacing
        spm_to_set = splrep(spm_geo, t_set_geo)
        t_set_spm_vals_interp = splev(spm_desired, spm_to_set)

        # ring longitude at final spacing
        spm_to_rhodot = splrep(spm_geo, rho_dot_kms_geo)
        rho_dot_kms_vals_interp = splev(spm_desired, spm_to_rhodot)

        # FILLERS FOR RADIUS CORRECTION
        rho_corr_pole_km_vals = np.zeros(len(spm_desired))
        rho_corr_timing_km_vals = np.zeros(len(spm_desired))

        # Set attributes
        self.rho_km_vals = rho_km_desired
        self.t_oet_spm_vals = spm_desired
        self.p_norm_vals = p_norm_vals
        self.phase_rad_vals = phase_rad_vals

        self.B_rad_vals = B_rad_vals_interp
        self.D_km_vals = D_km_vals_interp
        self.F_km_vals = F_km_vals_interp
        self.f_sky_hz_vals = f_sky_hz_vals_interp
        self.phi_rad_vals = phi_ora_rad_vals_interp
        self.t_ret_spm_vals = t_ret_spm_vals_interp
        self.t_set_spm_vals = t_set_spm_vals_interp
        self.phi_rl_rad_vals = phi_rl_rad_vals_interp
        self.rho_dot_kms_vals = rho_dot_kms_vals_interp
        self.rho_corr_pole_km_vals = rho_corr_pole_km_vals
        self.rho_corr_timing_km_vals = rho_corr_timing_km_vals
        self.tau_vals = -np.sin(B_rad_vals_interp)*np.log(p_norm_vals)
        self.raw_tau_threshold_vals = np.interp(spm_desired,self.spm_thresh,self.tau_thresh)
        self.tau_threshold_vals = np.interp(spm_desired,self.spm_thresh,self.tau_thresh)






    def interp_IQ(self, rho_km_vals, IQ_c_vals, dr_km, profile_range):
        rho_km_desired = np.arange(profile_range[0], profile_range[1],
            step = dr_km)

        Ic = np.real(IQ_c_vals)
        Qc = np.imag(IQ_c_vals)

        rho_to_Ic = splrep(rho_km_vals, Ic)
        Ic_desired = splev(rho_km_desired, rho_to_Ic)


        rho_to_Qc = splrep(rho_km_vals, Qc)
        Qc_desired = splev(rho_km_desired, rho_to_Qc)

        IQc_desired = Ic_desired + 1j * Qc_desired

        return rho_km_desired, IQc_desired

    @classmethod
    def create_dlps(cls, rsr_inst, geo_inst, cal_inst, dr_km, dr_km_tol=0.01,
        verbose=False, write_file=False, profile_range=[65000., 150000.]):
        # check for ingress, egress, or chord
        prof_dir = geo_inst.rev_info['prof_dir']

        if prof_dir == '"INGRESS"':
            dlp_egr = None
            dlp_ing = cls(rsr_inst, geo_inst, cal_inst, dr_km,
                    dr_km_tol=dr_km_tol, verbose=verbose,
                    write_file=write_file, profile_range=profile_range)

        elif prof_dir == '"EGRESS"':
            dlp_ing = None
            dlp_egr = cls(rsr_inst, geo_inst, cal_inst, dr_km,
                    dr_km_tol=dr_km_tol, verbose=verbose,
                    write_file=write_file, profile_range=profile_range)

        elif prof_dir == '"BOTH"':
            # split ingress and egress
            geo_ing = copy.deepcopy(geo_inst)
            geo_egr = copy.deepcopy(geo_inst)
            ind = geo_inst.split_ind
            ind_spm1 = geo_inst.t_oet_spm_vals[ind]

            rsr_ing = copy.deepcopy(rsr_inst)
            rsr_egr = copy.deepcopy(rsr_inst)

            cal_ing = copy.deepcopy(cal_inst)
            cal_egr = copy.deepcopy(cal_inst)

            ind_spm2 = np.argwhere(rsr_inst.spm_vals == ind_spm1)[0][0]
            ind_spm3 = np.argwhere(cal_inst.t_oet_spm_vals == ind_spm1)[0][0]


            # split raw spm and IQ attributes
            rsr_ing.spm_vals = rsr_inst.spm_vals[:ind_spm2]
            rsr_ing.IQ_m = rsr_inst.IQ_m[:ind_spm2]

            rsr_egr.spm_vals = rsr_inst.spm_vals[ind_spm2:]
            rsr_egr.IQ_m = rsr_inst.IQ_m[ind_spm2:]

            # split calibration attributes
            cal_ing.IQ_c = cal_inst.IQ_c[:ind_spm2]
            cal_ing.p_free_vals = cal_inst.p_free_vals[:ind_spm3]
            cal_ing.f_sky_hz_vals = cal_inst.f_sky_hz_vals[:ind_spm3]
            cal_ing.t_oet_spm_vals = cal_inst.t_oet_spm_vals[:ind_spm3]
            cal_ing.f_sky_resid_fit_vals = cal_inst.f_sky_resid_fit_vals[:ind_spm3]

            cal_egr.IQ_c = cal_inst.IQ_c[ind_spm2:]
            cal_egr.p_free_vals = cal_inst.p_free_vals[ind_spm3:]
            cal_egr.f_sky_hz_vals = cal_inst.f_sky_hz_vals[ind_spm3:]
            cal_egr.t_oet_spm_vals = cal_inst.t_oet_spm_vals[ind_spm3:]
            cal_egr.f_sky_resid_fit_vals = cal_inst.f_sky_resid_fit_vals[ind_spm3:]

            # split geometry attributes
            geo_ing.rho_km_vals = geo_inst.rho_km_vals[:ind]
            geo_ing.rho_dot_kms_vals = geo_inst.rho_dot_kms_vals[:ind]
            geo_ing.t_oet_spm_vals = geo_inst.t_oet_spm_vals[:ind]
            geo_ing.t_ret_spm_vals = geo_inst.t_ret_spm_vals[:ind]
            geo_ing.t_set_spm_vals = geo_inst.t_set_spm_vals[:ind]
            geo_ing.B_deg_vals = geo_inst.B_deg_vals[:ind]
            geo_ing.phi_ora_deg_vals = geo_inst.phi_ora_deg_vals[:ind]
            geo_ing.phi_rl_deg_vals = geo_inst.phi_rl_deg_vals[:ind]
            geo_ing.D_km_vals = geo_inst.D_km_vals[:ind]
            geo_ing.F_km_vals = geo_inst.F_km_vals[:ind]
            geo_ing.rev_info['prof_dir'] = '"INGRESS"'

            dlp_ing = cls(rsr_ing, geo_ing, cal_ing, dr_km,
                    dr_km_tol=dr_km_tol, verbose=verbose,
                    write_file=write_file, profile_range=profile_range)


            geo_egr.rho_km_vals = geo_inst.rho_km_vals[ind:]
            geo_egr.rho_dot_kms_vals = geo_inst.rho_dot_kms_vals[ind:]
            geo_egr.t_oet_spm_vals = geo_inst.t_oet_spm_vals[ind:]
            geo_egr.t_ret_spm_vals = geo_inst.t_ret_spm_vals[ind:]
            geo_egr.t_set_spm_vals = geo_inst.t_set_spm_vals[ind:]
            geo_egr.B_deg_vals = geo_inst.B_deg_vals[ind:]
            geo_egr.phi_ora_deg_vals = geo_inst.phi_ora_deg_vals[ind:]
            geo_egr.phi_rl_deg_vals = geo_inst.phi_rl_deg_vals[ind:]
            geo_egr.D_km_vals = geo_inst.D_km_vals[ind:]
            geo_egr.F_km_vals = geo_inst.F_km_vals[ind:]
            geo_ing.rev_info['prof_dir'] = '"EGRESS"'

            dlp_egr = cls(rsr_egr, geo_egr, cal_egr, dr_km,
                    dr_km_tol=dr_km_tol, verbose=verbose, write_file=write_file,
                    profile_range=profile_range)
        else:
            print('WARNING (create_dlps()): Invalid profile direction given!')

        return dlp_ing, dlp_egr

"""
Revisions:
"""
