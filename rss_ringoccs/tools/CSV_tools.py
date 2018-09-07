import numpy as np
import pandas as pd
from scipy import interpolate
from .write_history_dict import write_history_dict

def get_geo(geodata,verbose=True):
    if verbose:
        print("\tExtracting Geo Data...")
    dfg = pd.read_csv(geodata, delimiter=',',
        names=[
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
        )
    if verbose:
        print("\tGeo Data Complete.")
    return dfg

def get_cal(caldata,verbose=True):
    if verbose:
        print("\tExtracting Cal Data...")
    dfc = pd.read_csv(caldata, delimiter=',',
        names=[
            "spm_vals",
            "f_sky_pred_vals",
            "f_sky_resid_fit_vals",
            "p_free_vals"
            ]
        )
    if verbose:
        print("\tCal Data Complete.")
    return dfc

def get_dlp(dlpdata,verbose=True):
    if verbose:
        print("\tExtracting DLP Data...")
    dfd = pd.read_csv(
        dlpdata, delimiter=',',
        names=[
            "rho_km_vals",
            "rho_corr_pole_km_vals",
            "rho_corr_timing_km_vals",
            "phi_rl_deg_vals",
            "phi_ora_deg_vals",
            "raw_tau_vals",
            "phase_deg_vals",
            "raw_tau_threshold_vals",
            "t_oet_spm_vals",
            "t_ret_spm_vals",
            "t_set_spm_vals",
            "B_deg_vals"
        ]
    )
    if verbose:
        print("\tDLP Data Complete")
    return dfd

def get_tau(taudata,verbose=True):
    if verbose:
        print("\tExtracting Tau Data...")
    dft = pd.read_csv(taudata, delimiter=',',
        names=[
            "rho_km_vals",
            "rho_km_pole_corr_vals",
            "rho_km_offsett_vals",
            "phi_rl_deg_vals",
            "phi_ora_deg_vals",
            "raw_tau_vals",
            "phase_deg_vals",
            "raw_tau_threshold_vals",
            "spm_vals",
            "t_ret_spm_vals",
            "t_set_spm_vals",
            "B_deg_vals"
        ]
    )
    if verbose:
        print("\tTau Data Complete")
    return dft


class ExtractCSVData(object):
    """
        Class:
            csv_extract_data
        Purpose:
            Read three csv files (Geo, Cal, and Tau) and return
            an instance containing all necessary attributes to run
            diffraction correction. This instance can be fed directly
            into the diffraction_correction class.
        Variables:
            Geo:    A string that contains the location of
                    the requested Geo file.
            Cal:    A string that contains the location of
                    the requested Cal file.
            Tau:    A string that contains the location of
                    the requested Tau file.
        Attributes:
            rho_km_vals:        Ring radius, in kilometers.
            phi_rad_vals:       Ring azimuth angle, in radians.
            p_norm_vals:        Normalized raw power.
            phase_rad_vals:     Difracted phase, in radians.
            B_rad_vals:         Ring opening angle, in radians.
            D_km_vals:          Distance from the spacecraft to the
                                ring intercept point, in kilometers.
            f_sky_hz_vals:      The frequency of the recieved
                                singals, in Hertz.
            rho_dot_kms_vals:   The time derivative of the ring
                                radius, in kilometers per second.
    """

    def __init__(self,geodata,caldata,dlpdata,taudata=None,verbose=True):

        if (not isinstance(geodata,str)):
            raise TypeError("geodata must be a string: '/path/to/geodata'")
        if (not isinstance(caldata,str)):
            raise TypeError("caldata must be a string: '/path/to/caldata'")
        if (not isinstance(dlpdata,str)):
            raise TypeError("dlpdata must be a string: '/path/to/dlpdata'")
        if (not isinstance(verbose,bool)):
            raise TypeError("verbose must be a boolean: True/False")

        # Save inputs as attributes.
        self.geodata = geodata
        self.caldata = caldata
        self.dlpdata = dlpdata
        self.taudata = taudata

        if verbose:
            print("Extracting Data from CSV Files:")

        # Extract GEO, CAL, and DLP data.
        geo_dat = get_geo(geodata,verbose=verbose)
        cal_dat = get_cal(caldata,verbose=verbose)
        dlp_dat = get_dlp(dlpdata,verbose=verbose)

        self.__retrieve_variables(geo_dat,cal_dat,dlp_dat,verbose)
        self.__compute_variables(verbose)
        self.__interpolate_variables(verbose)
        self.__del_attributes()

        if (not isinstance(taudata,type(None))):
            if (not isinstance(taudata,str)):
                raise TypeError("taudata must be a string: '/path/to/taudata'")
            else:
                tau_dat         = get_tau(taudata,verbose=verbose)
                rho_km_vals     = self.rho_km_vals
                tr              = tau_dat.rho_km_vals
                tt              = tau_dat.raw_tau_vals
                tpdeg           = tau_dat.phase_deg_vals
                tbdeg           = tau_dat.B_deg_vals
                rmin            = np.min(tr)
                rmax            = np.max(tr)
                rstart          = int(np.min((rho_km_vals-rmin>=0).nonzero()))
                rfin            = int(np.max((rmax-rho_km_vals>=0).nonzero()))
                tau_rho         = rho_km_vals[rstart:rfin+1]
                tbrad           = np.deg2rad(tbdeg)
                tprad           = np.deg2rad(tpdeg)
                tm              = np.sin(np.abs(tbrad))
                tt_interp       = interpolate.interp1d(tr,tt,kind='linear')
                tphase_interp   = interpolate.interp1d(tr,tprad,kind='linear')
                phase_vals      = tphase_interp(tau_rho)
                tau_vals        = tt_interp(tau_rho)
                tm_interp       = interpolate.interp1d(tr,tm,kind='linear')
                tmu             = tm_interp(tau_rho)
                power_vals      = np.exp(-tau_vals/tmu)
                self.power_vals = power_vals
                self.tau_vals   = tau_vals
                self.phase_vals = phase_vals
                self.tau_rho    = tau_rho

        if verbose:
            print("\tData Extraction Complete.")
        if verbose:
            print("\tWriting History...")

        input_vars = {
            "GEO Data":self.geodata,
            "CAL Data":self.caldata,
            "DLP Data":self.dlpdata
            }
        input_kwds = {
            "TAU Data": self.taudata,
            "Use of Verbose": verbose
            }
        self.history = write_history_dict(input_vars, input_kwds, __file__)
        if verbose:
            print("\tHistory Complete.")
        if verbose:
            print("\tExtract CSV Data Complete.")

    def __retrieve_variables(self,geo_dat,cal_dat,dlp_dat,verbose):
        if verbose:
            print("\tRetrieving Variables...")

        # Run an error check on rho_km_vals
        rho_km_vals = np.array(dlp_dat.rho_km_vals)
        if not isinstance(rho_km_vals, np.ndarray):
            raise TypeError("Bad DLP: rho_km_vals must be a numpy array")
        elif (not np.isreal(rho_km_vals).all()):
            raise ValueError("Bad DLP: rho_km_vals must be real valued")
        elif (np.min(rho_km_vals) < 0.0):
            raise ValueError("Bad DLP: rho_km_vals has negative values")
        else:
            pass
        
        # Run an error check on phi_ora_deg_vals
        phi_ora_deg_vals = np.array(dlp_dat.phi_ora_deg_vals)
        if not isinstance(phi_ora_deg_vals, np.ndarray):
            raise TypeError("Bad DLP: phi_ora_deg_vals must be a numpy array")
        elif (not np.isreal(phi_ora_deg_vals).all()):
            raise ValueError("Bad DLP: phi_ora_deg_vals must be real valued")
        elif (np.max(np.abs(phi_ora_deg_vals)) > 360.0):
            raise ValueError("Bad DLP: max{|phi_ora_deg_vals|} > 360")
        else:
            pass

        # Run an error check on raw_tau_vals
        raw_tau_vals = np.array(dlp_dat.raw_tau_vals)
        if (not isinstance(raw_tau_vals, np.ndarray)):
            raise TypeError("Bad DLP: raw_tau_vals must be a numpy array")
        elif (not np.isreal(raw_tau_vals).all()):
            raise ValueError("raw_tau_vals must be real valued")
        else:
            pass

        # Run an error check on phase_deg_vals
        phase_deg_vals = np.array(dlp_dat.phase_deg_vals)
        if (not isinstance(phase_deg_vals,np.ndarray)):
            raise TypeError("Bad DLP: phase_deg_vals must be a numpy array")
        elif (not np.isreal(phase_deg_vals).all()):
            raise ValueError("Bad DLP: phase_deg_vals must be real valued")
        elif (np.max(np.abs(phase_deg_vals)) > 360.0):
            raise ValueError("Bad DLP: max{|phase_deg_vals|} > 360")
        else:
            pass
    
        # Run an error check on B_deg_vals
        B_deg_vals = np.array(dlp_dat.B_deg_vals)
        if (not isinstance(B_deg_vals,np.ndarray)):
            raise TypeError("Bad DLP: B_deg_vals must be a numpy array")
        elif (not np.isreal(B_deg_vals).all()):
            raise ValueError("Bad DLP: B_deg_vals must be real valued")
        elif (np.max(np.abs(phi_ora_deg_vals)) > 360.0):
            raise ValueError("Bad DLP: max{|B_deg_vals|} > 360")
        else:
            pass

        # Run an error check on t_ret_spm_vals
        t_ret_spm_vals = np.array(dlp_dat.t_ret_spm_vals)
        if (not isinstance(t_ret_spm_vals,np.ndarray)):
            raise TypeError("Bad DLP: t_ret_spm_vals must be a numpy array")
        elif (not np.isreal(t_ret_spm_vals).all()):
            raise ValueError("Bad DLP: t_ret_spm_vals must be real valued")
        elif (np.min(t_ret_spm_vals) < 0.0):
            raise ValueError("Bad DLP: t_ret_spm_vals has negative values.")
        else:
            pass

        # Run an error check on t_set_spm_vals
        t_set_spm_vals = np.array(dlp_dat.t_set_spm_vals)
        if (not isinstance(t_set_spm_vals,np.ndarray)):
            raise TypeError("Bad DLP: t_set_spm_vals must be a numpy array")
        elif (not np.isreal(t_set_spm_vals).all()):
            raise ValueError("Bad DLP: t_set_spm_vals must be real valued")
        elif (np.min(t_ret_spm_vals) < 0.0):
            raise ValueError("Bad DLP: t_set_spm_vals has negative values.")
        else:
            pass

        # Run an error check on t_oet_spm_vals
        t_oet_spm_vals = np.array(dlp_dat.t_oet_spm_vals)
        if (not isinstance(t_oet_spm_vals, np.ndarray)):
            raise TypeError("Bad DLP: t_oet_spm_vals must be a numpy array")
        elif (not np.isreal(t_oet_spm_vals).all()):
            raise ValueError("Bad DLP: t_oet_spm_vals must be real valued")
        elif (np.min(t_ret_spm_vals) < 0.0):
            raise ValueError("Bad DLP: t_oet_spm_vals has negative values")
        else:
            pass

        # Run an error check on rho_corr_pole_km_vals
        rho_corr_pole_km_vals = np.array(dlp_dat.rho_corr_pole_km_vals)
        if (not isinstance(rho_corr_pole_km_vals, np.ndarray)):
            raise TypeError(
                "Bad DLP: rho_corr_pole_km_vals must be a numpy array"
            )
        elif (not np.isreal(rho_corr_pole_km_vals).all()):
            raise ValueError(
                "Bad DLP: rho_corr_pole_km_vals must be real valued"
            )
        else:
            pass

        # Run an error check on rho_corr_timing_km_vals
        rho_corr_timing_km_vals = np.array(dlp_dat.rho_corr_timing_km_vals)
        if (not isinstance(rho_corr_timing_km_vals, np.ndarray)):
            raise TypeError(
                "Bad DLP: rho_corr_timing_km_vals must be a numpy array"
            )
        elif (not np.isreal(rho_corr_timing_km_vals).all()):
            raise ValueError(
                "Bad DLP: rho_corr_timing_km_vals must be real valued"
            )
        else:
            pass

        # Run an error check on phi_rl_deg_vals
        phi_rl_deg_vals = np.array(dlp_dat.phi_rl_deg_vals)
        if (not isinstance(phi_rl_deg_vals, np.ndarray)):
            raise TypeError("Bad DLP: phi_rl_deg_vals must be a numpy array")
        elif (not np.isreal(phi_rl_deg_vals).all()):
            raise ValueError("Bad DLP: phi_rl_deg_vals must be real valued")
        else:
            pass

        # Run an error check on raw_tau_threshold_vals
        raw_tau_threshold_vals = np.array(dlp_dat.raw_tau_threshold_vals)
        if (not isinstance(raw_tau_threshold_vals, np.ndarray)):
            raise TypeError(
                "Bad DLP: raw_tau_threshold_vals must be a numpy array"
            )
        elif (not np.isreal(raw_tau_threshold_vals).all()):
            raise ValueError(
                "Bad DLP: raw_tau_threshold_vals must be real valued"
            )
        else:
            pass

        # Run an error check on geo_rho
        geo_rho = np.array(geo_dat.rho_km_vals)
        if (not isinstance(geo_rho, np.ndarray)):
            raise TypeError("Bad GEO: rho_km_vals must be a numpy array")
        elif (not np.isreal(geo_rho).all()):
            raise ValueError("Bad GEO: rho_km_vals must be real valued")
        elif (np.min(geo_rho) < 0.0):
            raise ValueError("Bad GEO: rho_km_vals has negative values.")
        else:
            pass

        # Run an error check on geo_D
        geo_D   = np.array(geo_dat.D_km_vals)
        if (not isinstance(geo_D, np.ndarray)):
            raise TypeError("Bad GEO: D_km_vals must be a numpy array")
        elif (not np.isreal(geo_D).all()):
            raise ValueError("Bad GEO: D_km_vals must be real valued")
        elif (np.min(geo_D) < 0.0):
            raise ValueError("Bad GEO: D_km_vals has negative values.")
        else:
            pass

        # Run an error check on geo_drho
        geo_drho    = np.array(geo_dat.rho_dot_kms_vals)
        if (not isinstance(geo_drho, np.ndarray)):
            raise TypeError("Bad GEO: rho_dot_kms_vals must be a numpy array")
        elif (not np.isreal(geo_drho).all()):
            raise ValueError("Bad GEO: rho_dot_kms_vals must be real valued")
        else:
            pass

        # Run an error check on f_sky_raw_vals
        f_sky_raw_vals = np.array(cal_dat.f_sky_pred_vals)
        if (not isinstance(f_sky_raw_vals, np.ndarray)):
            raise TypeError("Bad CAL: f_sky_raw_vals must be a numpy array")
        elif (not np.isreal(f_sky_raw_vals).all()):
            raise ValueError("Bad CAL: f_sky_raw_vals must be real valued")
        elif (np.min(f_sky_raw_vals) < 0.0):
            raise ValueError("Bad CAL: f_sky_raw_vals has negative values.")
        elif (np.min(f_sky_raw_vals) < 10000.0):
            raise ValueError("Bad CAL: f_sky_raw_vals less than 1e4 Hz.")
        else:
            pass

        self.rho_km_vals             = rho_km_vals
        self.phi_ora_deg_vals        = phi_ora_deg_vals
        self.raw_tau_vals            = raw_tau_vals
        self.phase_deg_vals          = phase_deg_vals
        self.B_deg_vals              = B_deg_vals
        self.geo_rho                 = geo_rho
        self.geo_D                   = geo_D
        self.geo_drho                = geo_drho

        # Set time and frequency attributes.
        self.f_sky_raw_vals = f_sky_raw_vals
        self.t_oet_spm_vals = t_oet_spm_vals
        self.t_ret_spm_vals = t_ret_spm_vals
        self.t_set_spm_vals = t_set_spm_vals
        self.rho_corr_pole_km_vals   = rho_corr_pole_km_vals
        self.rho_corr_timing_km_vals = rho_corr_timing_km_vals
        self.phi_rl_deg_vals         = phi_rl_deg_vals
        self.raw_tau_threshold_vals  = raw_tau_threshold_vals

    def __compute_variables(self,verbose):
        if verbose:
            print("\tComputing Variables...")
        phi_rad_vals    = np.deg2rad(self.phi_ora_deg_vals)
        phi_rl_rad_vals = np.deg2rad(self.phi_rl_deg_vals)
        phase_rad_vals  = np.deg2rad(self.phase_deg_vals)
        B_rad_vals      = np.deg2rad(self.B_deg_vals)
        raw_mu          = np.sin(np.abs(B_rad_vals))
        raw_tau_vals    = self.raw_tau_vals
        geo_drho        = self.geo_drho
        p_norm_vals     = np.exp(-raw_tau_vals/raw_mu)

        rho_km_vals = self.rho_km_vals
        t_oet_spm_vals = self.t_oet_spm_vals

        if (np.size(rho_km_vals) != np.size(t_oet_spm_vals)):
            raise ValueError("len(rho_km_vals) != len(t_oet_spm_vals")

        dr = np.zeros(np.size(rho_km_vals) - 1)
        dt = np.zeros(np.size(t_oet_spm_vals)-1)

        for i in range(np.size(rho_km_vals) - 1):
            dr[i] = rho_km_vals[i+1] - rho_km_vals[i]
            dt[i] = t_oet_spm_vals[i+1] - t_oet_spm_vals[i]
        
        drdt = dr/dt

        if (np.min(drdt) < 0.0) and (np.max (drdt) > 0.0):
            raise ValueError(
                "Bad DLP: drho/dt has positive and negative values."
            )
        elif (np.size((drdt == 0).nonzero()) != 0):
            raise ValueError("Bad DLP: drho/dt has zero valued elements.")
        elif (drdt < 0.0).all():
            occ = 'ingress'
        elif (drdt > 0.0).all():
            occ = 'egress'
        else:
            raise ValueError("Bad DLP: drho/dt has incompatible elements")

        if (occ == 'ingress'):  crange = (geo_drho < 0.0).nonzero()
        elif (occ == 'egress'): crange = (geo_drho > 0.0).nonzero()
        else:
            crange_e = (geo_drho > 0.0).nonzero()
            crange_i = (geo_drho < 0.0).nonzero()
            n_e      = np.size(crange_e)
            n_i      = np.size(crange_i)
            if (n_e != 0) and (n_i !=0):
                raise ValueError(
                    "rho_dot_kms_vals has positive and negative values.\
                    This is likely a chord occultation. Set occ='ingress'\
                    to examine the ingress portion, and occ='egress'\
                    for the egress porition.")
            elif (n_e == 0) and (n_i == 0):
                raise ValueError("rho_dot_kms_vals is either empty or zero.")
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else: raise TypeError("Bad Input: GEO DATA")
            del n_e, n_i, crange_e, crange_i

        if (np.size(crange) == 0):
            if (occ == 'ingress'):
                mes = "rho_dot_kms_vals is never negative."
            elif (occ == 'egress'):
                mes = "rho_dot_kms_vals is never positive."
            else: raise ValueError("Bad occ input: Set 'egress' or 'ingress'")
            raise ValueError("Bad occ Input: '%s': %s" % (occ,mes))
        
        self.occ             = occ
        self.crange          = crange
        self.p_norm_vals     = p_norm_vals
        self.B_rad_vals      = B_rad_vals
        self.phase_rad_vals  = phase_rad_vals
        self.phi_rad_vals    = phi_rad_vals
        self.phi_rl_rad_vals = phi_rl_rad_vals
        self.raw_mu          = raw_mu

    def __interpolate_variables(self,verbose):
        if verbose: print("\tInterpolating Data...")
        crange                  = self.crange
        rho_km_vals             = self.rho_km_vals
        t_set_spm_vals          = self.t_set_spm_vals
        t_ret_spm_vals          = self.t_ret_spm_vals
        t_oet_spm_vals          = self.t_oet_spm_vals
        rho_corr_pole_km_vals      = self.rho_corr_pole_km_vals
        rho_corr_timing_km_vals = self.rho_corr_timing_km_vals

        f_sky_raw_vals   = self.f_sky_raw_vals
        geo_rho          = self.geo_rho[crange]
        geo_drho         = self.geo_drho[crange]
        geo_D            = self.geo_D[crange]
        rmin             = np.min(geo_rho)
        rmax             = np.max(geo_rho)
        rstart           = int(np.min((rho_km_vals-rmin>=0.0).nonzero()))
        rfin             = int(np.max((rmax-rho_km_vals>=0.0).nonzero()))
        rho_km_vals      = rho_km_vals[rstart:rfin+1]
        phi_rad_vals     = self.phi_rad_vals[rstart:rfin+1]
        phase_rad_vals   = self.phase_rad_vals[rstart:rfin+1]
        B_rad_vals       = self.B_rad_vals[rstart:rfin+1]
        p_norm_vals      = self.p_norm_vals[rstart:rfin+1]
        d_km_interp      = interpolate.interp1d(geo_rho,geo_D,kind='linear')
        D_km_vals        = d_km_interp(rho_km_vals)
        rho_dot_interp   = interpolate.interp1d(geo_rho,geo_drho,kind='linear')
        rho_dot_kms_vals = rho_dot_interp(rho_km_vals)
        n_rho_vals       = np.size(rho_km_vals)
        n_f_vals         = np.size(f_sky_raw_vals)
        frange           = np.arange(n_f_vals)
        xrange           = np.arange(n_rho_vals)*(n_f_vals-1.0)/(n_rho_vals-1.0)
        fsky_interp      = interpolate.interp1d(frange,f_sky_raw_vals,kind='linear')
        f_sky_hz_vals    = fsky_interp(xrange)

        t_ret_spm_vals          = t_ret_spm_vals[rstart:rfin+1]
        t_set_spm_vals          = t_set_spm_vals[rstart:rfin+1]
        t_oet_spm_vals          = t_oet_spm_vals[rstart:rfin+1]
        rho_corr_pole_km_vals      = rho_corr_pole_km_vals[rstart:rfin+1]
        rho_corr_timing_km_vals = rho_corr_timing_km_vals[rstart:rfin+1]

        del f_sky_raw_vals,geo_rho,geo_drho,geo_D,rmin,rmax,rstart,rfin
        del rho_dot_interp,n_rho_vals,n_f_vals,frange,xrange,fsky_interp
        
        self.rho_km_vals        = rho_km_vals
        self.phi_rad_vals       = phi_rad_vals
        self.p_norm_vals        = p_norm_vals
        self.phase_rad_vals     = phase_rad_vals
        self.B_rad_vals         = B_rad_vals
        self.D_km_vals          = D_km_vals
        self.f_sky_hz_vals      = f_sky_hz_vals
        self.rho_dot_kms_vals   = rho_dot_kms_vals
        t_set_spm_vals          = t_set_spm_vals
        t_ret_spm_vals          = t_ret_spm_vals
        t_oet_spm_vals          = t_oet_spm_vals
        rho_corr_pole_km_vals   = rho_corr_pole_km_vals
        rho_corr_timing_km_vals = rho_corr_timing_km_vals

    def __del_attributes(self):
        del self.phi_ora_deg_vals,self.raw_tau_vals
        del self.phase_deg_vals,self.raw_mu,self.B_deg_vals
        del self.geo_rho,self.geo_D,self.geo_drho,self.crange


class GetUranusData(object):
    def __init__(self,geodata,dlpdata,dx=0.25,occ=None,verbose=False):
        if (not isinstance(geodata,str)):
            raise TypeError("geodata must be a string: '/path/to/geodata'")
        if (not isinstance(dlpdata,str)):
            raise TypeError("dlpdata must be a string: '/path/to/dlpdata'")
        if (not isinstance(dx,float)):
            raise TypeError("dx must be a floating point number")
        if (dx <= 0.0):
            raise ValueEorr("dx must be a positive number")
        if occ:
            if (not isinstance(occ,str)):
                raise TypeError("occ must be a string")
            else:
                occ = occ.replace(" ", "").lower()
                if (occ != 'ingress') and (occ != 'egress'):
                    raise ValueError("occ must be 'egress' of 'ingress'")
                else:
                    pass

        geo_dat = get_geo(geodata,verbose=verbose)
        dlp_dat = pd.read_csv(
            dlpdata, delimiter=',',
            names=[
                "t_oet_spm_vals",
                "p_norm_vals",
                "phase_rad_vals",
                "f_sky_hz_vals"
            ]
        )

        dlp_spm = np.array(dlp_dat.t_oet_spm_vals)
        dlp_pow = np.array(dlp_dat.p_norm_vals)
        dlp_phs = np.array(dlp_dat.phase_rad_vals)
        dlp_frq = np.array(dlp_dat.f_sky_hz_vals)
        geo_spm = np.array(geo_dat.t_oet_spm_vals)

        geo_rho = geo_dat.rho_km_vals
        n_rho = np.size(geo_rho)
        drho = np.zeros(n_rho-1)
        geo_D = geo_dat.D_km_vals
        geo_B = geo_dat.B_deg_vals
        geo_drho = geo_dat.rho_dot_kms_vals
        geo_phi = geo_dat.phi_ora_deg_vals
        t_dlp1 = np.min(dlp_spm)
        t_dlp2 = np.max(dlp_spm)
        t_geo1 = np.min(geo_spm)
        t_geo2 = np.max(geo_spm)

        t1 = np.max([t_dlp1,t_geo1])
        t2 = np.min([t_dlp2,t_geo2])
        if (t1 > t2):
            raise ValueError(
                "Geo and DLP data never overlap. No data available."
            )

        start = np.min((geo_spm >= t1).nonzero())
        finish = np.max((geo_spm <= t2).nonzero())
        tstart = np.min((dlp_spm >= t1).nonzero())
        tfinish = np.max((dlp_spm <= t2).nonzero())
        t_dlp = dlp_spm[tstart:tfinish+1]
        
        for i in range(n_rho-1):
            drho[i] = geo_rho[i+1]-geo_rho[i]
        if not occ:
            if (np.min(drho) < 0.0) and (np.max(drho) > 0.0):
                raise ValueError(
                    "\n\tdrho is positive and negative.\n\
                     \tSet occ to ingress or egress"
                )
            elif (drho > 0).all():
                crange = np.arange(start,finish+1)
            elif (drho < 0).all():
                crange = np.arange(start,finish+1)
                crange = crange[::-1]
            elif (drho == 0).all():
                raise ValueError("drho/dt = 0 for all points.")
            else:
                raise ValueError("drho/dt has invalid values.")
        elif (occ == 'ingress'):
            crange = (drho < 0.0).nonzero()
            if (np.size(crange) == 0):
                raise TypeError("drho is never negative. Use occ = 'egress'")
        elif (occ == 'egress'):
            crange = (drho > 0.0).nonzero()
            if (np.size(crange) == 0):
                raise TypeError("drho is never positive. Use occ = 'ingress'")
        else:
            raise ValueError("Invalid occ keyword: %s" % occ)

        dlp_rho_interp = interpolate.interp1d(geo_spm,geo_rho,kind="linear")
        dlp_rho = dlp_rho_interp(t_dlp)

        rho_min = np.min(dlp_rho)
        rho_max = np.max(dlp_rho)
        self.rho_km_vals = np.arange(rho_min,rho_max,dx)

        drho_interp = interpolate.interp1d(geo_rho,geo_drho,kind="linear")
        self.rho_dot_kms_vals = drho_interp(self.rho_km_vals)
        
        geo_D_interp = interpolate.interp1d(geo_rho,geo_D,kind="linear")
        self.D_km_vals = geo_D_interp(self.rho_km_vals)

        geo_B_interp = interpolate.interp1d(geo_rho,geo_B,kind="linear")
        self.B_rad_vals = np.deg2rad(geo_B_interp(self.rho_km_vals))

        geo_phi_interp = interpolate.interp1d(geo_rho,geo_phi,kind="linear")
        self.phi_rad_vals = np.deg2rad(geo_phi_interp(self.rho_km_vals))

        power_interp = interpolate.interp1d(dlp_rho,dlp_pow,kind="linear")
        self.p_norm_vals = power_interp(self.rho_km_vals)

        phase_interp = interpolate.interp1d(dlp_rho,dlp_phs[tstart:tfinish+1],kind="linear")
        self.phase_rad_vals = phase_interp(self.rho_km_vals)

        freq_interp = interpolate.interp1d(dlp_rho,dlp_frq[tstart:tfinish+1],kind="linear")
        self.f_sky_hz_vals = freq_interp(self.rho_km_vals)

        n = np.size(self.rho_km_vals)

        self.t_oet_spm_vals = np.zeros(n)
        self.t_ret_spm_vals = np.zeros(n)
        self.t_set_spm_vals = np.zeros(n)
        self.rho_corr_pole_km_vals = np.zeros(n)
        self.rho_corr_timing_km_vals = np.zeros(n)
        self.phi_rl_rad_vals = np.zeros(n)
        self.raw_tau_threshold_vals = np.zeros(n)

        if verbose:
            print("\tData Extraction Complete.")
        if verbose:
            print("\tWriting History...")

        input_vars = {
            "GEO Data": geodata,
            "DLP Data": dlpdata
            }
        input_kwds = {
            "occ":             occ,
            "Use of Verbose":  verbose
            }
        self.history = write_history_dict(input_vars, input_kwds, __file__)
        if verbose:
            print("\tHistory Complete.")
        if verbose:
            print("\tExtract CSV Data Complete.")


class PureCSVReader(object):
    def __init__(self,dat):
        if (not isinstance(dat, str)):
            raise TypeError("Text file must be a string.")
        df = pd.read_csv(dat)
        self.rho_km_vals      = np.array(df.rho_km_vals)
        self.phase_rad_vals   = np.array(df.phase_rad_vals)
        self.p_norm_vals      = np.array(df.p_norm_vals)
        self.phi_rad_vals     = np.array(df.phi_rad_vals)
        self.B_rad_vals       = np.array(df.B_rad_vals)
        self.f_sky_hz_vals    = np.array(df.f_sky_hz_vals)
        self.D_km_vals        = np.array(df.D_km_vals)
        self.rho_dot_kms_vals = np.array(df.rho_dot_kms_vals)
