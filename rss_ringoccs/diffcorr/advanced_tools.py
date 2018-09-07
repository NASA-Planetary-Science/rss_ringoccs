import sys
import numpy as np
from scipy import interpolate
from .diffraction_correction import DiffractionCorrection
from .physics_functions import fresnel_forward, freq_wav
from rss_ringoccs.tools.misc_tools import get_geo, ExtractCSVData

class compare_tau(object):
    def __init__(self,geo,cal,dlp,tau,res,rng='all',
        wtype="kb25",bfac=True,fft=False,verbose=True,norm=True):

        data        = ExtractCSVData(
            geo,cal,dlp,taudata=tau,verbose=verbose
        )
        tr          = data.tau_rho
        tau_power   = data.power_vals
        tau_phase   = data.phase_vals
        tau_tau     = data.tau_vals
        rec         = DiffractionCorrection(
            data,res,rng=rng,bfac=bfac,wtype=wtype,fft=fft,
            verbose=verbose,norm=norm
        )
        rho_km_vals = rec.rho_km_vals
        rmin,rmax = np.min(rho_km_vals),np.max(rho_km_vals)
        tau_rstart  = int(np.min((tr-rmin>=0).nonzero()))
        tau_rfin    = int(np.max((rmax-tr>=0).nonzero()))

        self.rho_km_vals = rec.rho_km_vals
        self.power_vals  = rec.power_vals
        self.tau_vals    = rec.tau_vals
        self.phase_vals  = rec.phase_vals
        self.tau_power   = tau_power[tau_rstart:tau_rfin+1]
        self.tau_tau     = tau_tau[tau_rstart:tau_rfin+1]
        self.tau_phase   = tau_phase[tau_rstart:tau_rfin+1]
        self.wtype       = wtype
        self.res         = res

class find_optimal_resolution(object):
    def __init__(self,geo,cal,dlp,tau,sres,dres,nres,
        rng='all',wlst=['kb25'],verbose=True):
        self.linfint = None
        self.l2int   = None
        self.linffft = None
        self.l2fft   = None

        nwins   = np.size(wlst)
        linfint = np.zeros((nwins,nres))
        l2int   = np.zeros((nwins,nres))
        linffft = np.zeros((nwins,nres))
        l2fft   = np.zeros((nwins,nres))
        resint  = np.zeros((nwins))
        resfft  = np.zeros((nwins))
        eres    = sres + (nres-1)*dres

        res         = sres
        data        = extract_csv_data(geo,cal,dlp,taudata=tau,verbose=verbose)
        tr          = data.tau_rho
        tau_power   = data.power_vals
        rec         = DiffractionCorrection(data,res,rng=rng,wtype="kb35",verbose=False)
        rho_km_vals = rec.rho_km_vals
        power_vals  = rec.power_vals
        rmin        = np.min(rho_km_vals)
        rmax        = np.max(rho_km_vals)
        rstart      = int(np.min((tr-rmin>=0).nonzero()))
        rfin        = int(np.max((rmax-tr>=0).nonzero()))
        tau_power   = tau_power[rstart:rfin+1]
        sys.stdout.write("\033[K")
        for i in np.arange(nres):
            for j in range(nwins):
                wtype        = wlst[j]       
                recint       = DiffractionCorrection(data,res,rng=rng,wtype=wtype,fft=False)
                p_int        = recint.power_vals
                linf         = np.max(np.abs(p_int - tau_power))
                l2           = np.sqrt(np.sum(np.abs(p_int-tau_power)**2)*recint.dx_km)
                linfint[j,i] = linf
                l2int[j,i]   = l2

                recfft       = DiffractionCorrection(data,res,rng=rng,wtype=wtype,fft=True)
                p_fft        = recfft.power_vals
                linf         = np.max(np.abs(p_fft - tau_power))
                l2           = np.sqrt(np.sum(np.abs(p_fft-tau_power)**2)*recint.dx_km)
                linffft[j,i] = linf
                l2fft[j,i]   = l2
                sys.stdout.write("\033[K")
                if verbose:
                    printmes = ('Res:',res,'Max:',eres,"WTYPE:",wtype)
                    print("%s %f %s %f %s %s" % printmes)
            res += dres
        for j in range(nwins):
            resint[j] = sres+dres*np.min((linfint[j,...] == np.min(linfint[j,...])).nonzero())
            resfft[j] = sres+dres*np.min((linffft[j,...] == np.min(linffft[j,...])).nonzero())
        self.linfint = linfint
        self.l2int   = l2int
        self.linffft = linffft
        self.l2fft   = l2fft
        self.resint  = resint
        self.resfft  = resfft

class delta_impulse_diffraction(object):
    def __init__(self,geo,lambda_km,res,rho,dx_km_desired=0.25,
        occ=False,wtype='kb25',fwd=False,norm=True,bfac=True,
        verbose=True,fft=False,psitype='full',usefres=False):

        if (not isinstance(res, float)) or (not isinstance(res, int)):
            try:
                res = float(res)
            except TypeError:
                raise TypeError("res must be a positive floating point number")
        if (res <= 0.0):
            raise ValueError("res must be a positive floating point number")
        if not isinstance(wtype, str):
            raise TypeError("wtype must be a string. Ex: 'coss'")
        if not isinstance(fwd, bool):
            raise TypeError("fwd must be Boolean: True/False")
        if not isinstance(norm, bool):
            raise TypeError("norm must be Boolean: True/False")
        if not isinstance(bfac, bool):
            raise TypeError("bfac must be Boolean: True/False")
        if not isinstance(fft, bool):
            raise TypeError("fft must be Boolean: True/False")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be Boolean: True/False")
        if not isinstance(psitype, str):
            raise TypeError("psitype must be a string. Ex: 'full'")
        if not isinstance(verbose, bool):
            raise TypeError("usefres must be Boolean: True/False")
    
        data = get_geo(geo,verbose=verbose)
        self.__retrieve_variables(data,verbose)
        self.__compute_variables(dx_km_desired,occ,verbose)
        self.__interpolate_variables(verbose)

        r0       = self.rho_km_vals
        b        = self.B_rad_vals
        d        = self.D_km_vals
        phi      = self.phi_rad_vals
        rng      = [rho-10.0*res,rho+10.0*res]
        nstar    = np.min((r0-rho>=0).nonzero())
        r        = r0[nstar]
        phistar  = phi[nstar]
        F        = fresnel_scale(lambda_km,d,phi,b)

        if (usefres == True):
            psi_vals = (np.pi/2.0)*((r0-r)/F)*((r0-r)/F)
        else:
            kD          = 2.0 * np.pi * d[nstar] / lambda_km
            dphi_s  = psi_factor(r0,r,b,phi)
            phi_s       = phi - dphi_s
            loop        = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s)) > 1.e-8):
                psi_d1  = psi_d1_phi(r,r0,d,b,phi_s,phi0,error_check=False)
                psi_d2  = psi_d2_phi(r,r0,d,b,phi_s,phi0,error_check=False)
                dphi_s  = -psi_d1 / psi_d2
                phi_s  += dphi_s
                loop   += 1
                if loop > 20:
                    break
                if verbose: print("Psi Iter: %d" % loop,end="\r")
            if verbose: print("Psi Iter: %d  dphi: %" % (loop,np.max(dphi_s)))
            psi_vals = kD * psi(r,r0,d,b,phi_s,phi0)

        T_hat_vals     = (1.0-1.0j)*np.exp(1j*psi_vals)/(2.0*F)
        p_norm_vals    = np.abs(T_hat_vals)*np.abs(T_hat_vals)
        phase_rad_vals = -np.arctan2(np.imag(T_hat_vals),np.real(T_hat_vals))

        lambda_vals         = np.zeros(np.size(self.rho_km_vals))+lambda_km
        self.p_norm_vals    = p_norm_vals
        self.phase_rad_vals = phase_rad_vals
        self.f_sky_hz_vals  = freq_wav(lambda_vals)
        self.nstar          = nstar
        self.history        = "Bob"

        recdata = diffraction_correction(self,res,rng=rng,
            wtype=wtype,fwd=fwd,norm=norm,verbose=verbose,
            bfac=bfac,fft=fft,psitype=psitype)

        self.rho_star           = rho
        self.rho_km_vals        = recdata.rho_km_vals
        self.p_norm_vals        = recdata.p_norm_vals
        self.phase_rad_vals     = recdata.phase_rad_vals
        self.B_rad_vals         = recdata.B_rad_vals
        self.D_km_vals          = recdata.D_km_vals
        self.f_sky_hz_vals      = recdata.f_sky_hz_vals
        self.phi_rad_vals       = recdata.phi_rad_vals
        self.rho_dot_kms_vals   = recdata.rho_dot_kms_vals
        self.T_hat_vals         = recdata.T_hat_vals
        self.res                = recdata.res
        self.wtype              = recdata.wtype
        self.rng                = recdata.rng
        self.F_km_vals          = recdata.F_km_vals
        self.w_km_vals          = recdata.w_km_vals
        self.mu_vals            = recdata.mu_vals
        self.lambda_sky_km_vals = recdata.lambda_sky_km_vals
        self.dx_km              = recdata.dx_km
        self.norm_eq            = recdata.norm_eq
        self.n_used             = recdata.n_used
        self.start              = recdata.start
        self.T_vals             = recdata.T_vals
        self.power_vals         = recdata.power_vals
        self.tau_vals           = recdata.tau_vals
        self.phase_vals         = recdata.phase_vals
        self.p_norm_fwd_vals    = recdata.p_norm_fwd_vals
        self.T_hat_fwd_vals     = recdata.T_hat_fwd_vals
        self.phase_fwd_vals     = recdata.phase_fwd_vals
        self.norm               = recdata.norm
        self.fwd                = recdata.fwd
        self.fft                = recdata.fft

    def __retrieve_variables(self,geo_dat,verbose):
        if verbose: print("Retrieving Variables...")
        geo_rho          = np.array(geo_dat.rho_km_vals)
        geo_D            = np.array(geo_dat.D_km_vals)
        geo_drho         = np.array(geo_dat.rho_dot_kms_vals)
        geo_phi          = np.array(geo_dat.phi_ora_deg_vals)
        geo_B            = np.array(geo_dat.B_deg_vals)
        rhotype,drhotype = check_real(geo_rho),check_real(geo_drho)
        phitype,Btype    = check_real(geo_phi),check_real(geo_B)
        Dtype            = check_real(geo_D)
        if not rhotype:
            raise TypeError("Bad GEO: rho_km_vals not real valued.")
        elif (np.min(geo_rho) < 0.0):
            raise ValueError("Bad GEO: rho_km_vals has negative values.")
        else: self.geo_rho  = geo_rho
        if not Dtype:
            raise TypeError("Bad GEO: D_km_vals not real valued.")
        elif (np.min(geo_D) < 0.0):
            raise ValueError("Bad GEO: D_km_vals has negative values.")
        else: self.geo_D    = geo_D
        if not drhotype:
            raise TypeError("Bad GEO: rho_dot_kms_vals not real valued.")
        else: self.geo_drho = geo_drho
        if not phitype:
            raise TypeError("Bad GEO: phi_deg_ora_vals not real valued.")
        elif (np.max(np.abs(geo_phi)) > 360.0):
            raise ValueError("Bad GEO: phi_deg_ora_vals > 360")
        else: self.geo_phi  = geo_phi
        if not Btype:
            raise TypeError("Bad GEO: B_deg_vals not real valued.")
        elif (np.max(np.abs(geo_B)) > 360.0):
            raise ValueError("Bad GEO: B_de2g_vals > 360")
        else: self.geo_B    = geo_B
        del rhotype,Dtype,drhotype,phitype,Btype

    def __compute_variables(self,dx_km_desired,occ,verbose):
        if verbose: print("Computing Variables...")
        geo_drho        = self.geo_drho

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
                    for the egress portion.")
            elif (n_e == 0) and (n_i == 0):
                raise ValueError("rho_dot_kms_vals is either empty or zero.")
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else: raise TypeError("Bad Input: GEO DATA")
        del n_e,n_i,crange_e,crange_i
        if (np.size(crange) == 0):
            if (occ == 'ingress'):
                mes = "rho_dot_kms_vals is never negative."
            elif (occ == 'egress'):
                mes = "rho_dot_kms_vals is never positive."
            else: raise ValueError("Bad occ input: Set 'egress' or 'ingress'")
            raise ValueError("Bad occ Input: '%s': %s" % (occ,mes))
        self.geo_rho        = self.geo_rho[crange]
        self.geo_D          = self.geo_D[crange]
        self.geo_drho       = self.geo_drho[crange]
        self.geo_phi        = self.geo_phi[crange]
        self.geo_B          = self.geo_B[crange]
        rmin                = np.min(geo_rho)
        rmax                = np.max(geo_rho)
        self.rho_km_vals    = np.arange(rmin,rmax,dx_km_desired)
        del rmin,rmax,geo_rho,geo_D,geo_drho,geo_phi,geo_B,rho_km_vals

    def __interpolate_variables(self,verbose):
        if verbose: print("Interpolating Data...")
        rho_km_vals           = self.rho_km_vals
        geo_rho               = self.geo_rho
        geo_drho              = self.geo_drho
        geo_D                 = self.geo_D
        geo_phi               = self.geo_phi
        geo_B                 = self.geo_B
        D_km_interp           = interpolate.interp1d(geo_rho,geo_D)
        self.D_km_vals        = D_km_interp(rho_km_vals)
        rho_dot_interp        = interpolate.interp1d(geo_rho,geo_drho)
        self.rho_dot_kms_vals = rho_dot_interp(rho_km_vals)
        phi_deg_interp        = interpolate.interp1d(geo_rho,geo_phi)
        phi_deg_vals          = phi_deg_interp(rho_km_vals)
        B_deg_interp          = interpolate.interp1d(geo_rho,geo_B)
        B_deg_vals            = B_deg_interp(rho_km_vals)
        self.phi_rad_vals     = np.deg2rad(phi_deg_vals)
        self.B_rad_vals       = np.deg2rad(B_deg_vals)
        del geo_rho,geo_drho,geo_D,geo_phi,geo_B,D_km_interp,rho_dot_interp
        del phi_deg_vals,phi_deg_interp,B_deg_interp,B_deg_vals

class double_slit_diffraction(object):
    def __init__(
        self,geo,lambda_km,res,rho_s,rho_e,rho1,rho2,rho3,rho4,dx_km_desired=0.25,
        occ=False,wtype='kb25',fwd=False,norm=True,bfac=True,
        verbose=True,fft=False,psitype='full',usefres=False):

        if (not isinstance(res, float)) or (not isinstance(res, int)):
            try:
                res = float(res)
            except TypeError:
                raise TypeError("res must be a positive floating point number")
        if (res <= 0.0):
            raise ValueError("res must be a positive floating point number")
        if not isinstance(wtype, str):
            raise TypeError("wtype must be a string. Ex: 'coss'")
        if not isinstance(fwd, bool):
            raise TypeError("fwd must be Boolean: True/False")
        if not isinstance(norm, bool):
            raise TypeError("norm must be Boolean: True/False")
        if not isinstance(bfac, bool):
            raise TypeError("bfac must be Boolean: True/False")
        if not isinstance(fft, bool):
            raise TypeError("fft must be Boolean: True/False")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be Boolean: True/False")
        if not isinstance(psitype, str):
            raise TypeError("psitype must be a string. Ex: 'full'")
        if not isinstance(verbose, bool):
            raise TypeError("usefres must be Boolean: True/False")
    
        data = get_geo(geo,verbose=verbose)
        self.__retrieve_variables(data,verbose)
        self.__compute_variables(dx_km_desired,occ,verbose)
        self.__interpolate_variables(verbose)

        r0       = self.rho_km_vals
        B        = self.B_rad_vals
        D        = self.D_km_vals
        phi      = self.phi_rad_vals
        rng      = [rho1-10.0*res,rho4+10.0*res]
        cb = np.cos(self.B_rad_vals)
        sb = np.sin(self.B_rad_vals)
        sp = np.sin(self.phi_rad_vals)

        F = np.sqrt(0.5 * lambda_km * D * (1 - (cb*cb) * (sp*sp)) / (sb*sb))
        w = 2*F*F/res
        w_max = np.max(w)
        nw = w_max/(2.0*dx_km_desired)

        T_vals = np.zeros(np.size(r0)) + 0.0j

        n1 = np.min((r0>=rho1).nonzero())
        n2 = np.max((r0<=rho2).nonzero())
        n3 = np.min((r0>=rho3).nonzero())
        n4 = np.max((r0<=rho4).nonzero())
        n_start = np.min((r0>rho_s-nw/2.0).nonzero())
        n_fin = np.max((r0<rho_e+nw/2.0).nonzero())
        n_used = n_fin - n_start + 1

        T_vals[n1:n2] = 1.0
        T_vals[n3:n4] = 1.0

        T_hat_vals = fresnel_forward(
            r0,F,phi,B,D,T_vals,lambda_km,w,dx_km_desired,wtype,n_start,
            n_used,norm=norm,fft=fft,verbose=verbose,psitype=psitype
        )

        p_norm_vals    = np.abs(T_hat_vals)*np.abs(T_hat_vals)
        phase_rad_vals = -np.arctan2(np.imag(T_hat_vals),np.real(T_hat_vals))

        lambda_vals         = np.zeros(np.size(self.rho_km_vals))+lambda_km
        self.p_norm_vals    = p_norm_vals
        self.phase_rad_vals = phase_rad_vals
        self.f_sky_hz_vals  = freq_wav(lambda_vals)
        self.nstar          = nstar
        self.history        = "Bob"

        recdata = diffraction_correction(self,res,rng=rng,
            wtype=wtype,fwd=fwd,norm=norm,verbose=verbose,
            bfac=bfac,fft=fft,psitype=psitype)

        self.rho_star           = rho
        self.rho_km_vals        = recdata.rho_km_vals
        self.p_norm_vals        = recdata.p_norm_vals
        self.phase_rad_vals     = recdata.phase_rad_vals
        self.B_rad_vals         = recdata.B_rad_vals
        self.D_km_vals          = recdata.D_km_vals
        self.f_sky_hz_vals      = recdata.f_sky_hz_vals
        self.phi_rad_vals       = recdata.phi_rad_vals
        self.rho_dot_kms_vals   = recdata.rho_dot_kms_vals
        self.T_hat_vals         = recdata.T_hat_vals
        self.res                = recdata.res
        self.wtype              = recdata.wtype
        self.rng                = recdata.rng
        self.F_km_vals          = recdata.F_km_vals
        self.w_km_vals          = recdata.w_km_vals
        self.mu_vals            = recdata.mu_vals
        self.lambda_sky_km_vals = recdata.lambda_sky_km_vals
        self.dx_km              = recdata.dx_km
        self.norm_eq            = recdata.norm_eq
        self.n_used             = recdata.n_used
        self.start              = recdata.start
        self.T_vals             = recdata.T_vals
        self.power_vals         = recdata.power_vals
        self.tau_vals           = recdata.tau_vals
        self.phase_vals         = recdata.phase_vals
        self.p_norm_fwd_vals    = recdata.p_norm_fwd_vals
        self.T_hat_fwd_vals     = recdata.T_hat_fwd_vals
        self.phase_fwd_vals     = recdata.phase_fwd_vals
        self.norm               = recdata.norm
        self.fwd                = recdata.fwd
        self.fft                = recdata.fft

    def __retrieve_variables(self,geo_dat,verbose):
        if verbose: print("Retrieving Variables...")
        geo_rho          = np.array(geo_dat.rho_km_vals)
        geo_D            = np.array(geo_dat.D_km_vals)
        geo_drho         = np.array(geo_dat.rho_dot_kms_vals)
        geo_phi          = np.array(geo_dat.phi_ora_deg_vals)
        geo_B            = np.array(geo_dat.B_deg_vals)

        # Run an error check on geo_rho
        geo_rho = np.array(geo_dat.rho_km_vals)
        if (not isinstance(geo_rho, np.ndarray)):
            raise TypeError("Bad GEO: rho_km_vals must be a numpy array")
        elif (not np.isreal(geo_rho).all()):
            raise ValueError("Bad GEO: rho_km_vals must be real valued")
        elif (np.min(geo_rho) < 0.0):
            raise ValueError("Bad GEO: rho_km_vals has negative values.")
        else:
            self.geo_rho = geo_rho

        # Run an error check on geo_D
        geo_D   = np.array(geo_dat.D_km_vals)
        if (not isinstance(geo_D, np.ndarray)):
            raise TypeError("Bad GEO: D_km_vals must be a numpy array")
        elif (not np.isreal(geo_D).all()):
            raise ValueError("Bad GEO: D_km_vals must be real valued")
        elif (np.min(geo_D) < 0.0):
            raise ValueError("Bad GEO: D_km_vals has negative values.")
        else:
            self.geo_D = geo_D

        # Run an error check on geo_drho
        geo_drho    = np.array(geo_dat.rho_dot_kms_vals)
        if (not isinstance(geo_drho, np.ndarray)):
            raise TypeError("Bad GEO: rho_dot_kms_vals must be a numpy array")
        elif (not np.isreal(geo_drho).all()):
            raise ValueError("Bad GEO: rho_dot_kms_vals must be real valued")
        else:
            self.geo_drho = geo_drho

        if (not isinstance(geo_phi,np.ndarray)):
            raise TypeError("Bad GEO: phi_deg_ora_vals must be a numpy array")
        elif (not np.isreal(geo_phi).all()):
            raise ValueError("Bad GEO: phi_deg_ora_vals must be real valued")
        elif (np.max(np.abs(geo_phi)) > 360.0):
            raise ValueError("Bad GEO: max{|phi_deg_ora_vals|} > 360")
        else:
            self.geo_phi  = geo_phi

        if (not isinstance(geo_B,np.ndarray)):
            raise TypeError("Bad GEO: B_deg_vals must be a numpy array")
        elif (not np.isreal(geo_B).all()):
            raise ValueError("Bad GEO: B_deg_vals must be real valued")
        elif (np.max(np.abs(geo_B)) > 360.0):
            raise ValueError("Bad GEO: max{|B_deg_vals|} > 360")
        else:
            self.geo_B  = geo_B

    def __compute_variables(self,dx_km_desired,occ,verbose):
        if verbose: print("Computing Variables...")
        geo_drho        = self.geo_drho

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
                    for the egress portion.")
            elif (n_e == 0) and (n_i == 0):
                raise ValueError("rho_dot_kms_vals is either empty or zero.")
            elif (n_e != 0) and (n_i == 0):
                crange = crange_e
                occ    = 'egress'
            elif (n_e == 0) and (n_i != 0):
                crange = crange_i
                occ    = 'ingress'
            else: raise TypeError("Bad Input: GEO DATA")
        del n_e,n_i,crange_e,crange_i
        if (np.size(crange) == 0):
            if (occ == 'ingress'):
                mes = "rho_dot_kms_vals is never negative."
            elif (occ == 'egress'):
                mes = "rho_dot_kms_vals is never positive."
            else: raise ValueError("Bad occ input: Set 'egress' or 'ingress'")
            raise ValueError("Bad occ Input: '%s': %s" % (occ,mes))
        self.geo_rho        = self.geo_rho[crange]
        self.geo_D          = self.geo_D[crange]
        self.geo_drho       = self.geo_drho[crange]
        self.geo_phi        = self.geo_phi[crange]
        self.geo_B          = self.geo_B[crange]
        rmin                = np.min(self.geo_rho)
        rmax                = np.max(self.geo_rho)
        self.rho_km_vals    = np.arange(rmin,rmax,dx_km_desired)

    def __interpolate_variables(self,verbose):
        if verbose: print("Interpolating Data...")
        rho_km_vals           = self.rho_km_vals
        geo_rho               = self.geo_rho
        geo_drho              = self.geo_drho
        geo_D                 = self.geo_D
        geo_phi               = self.geo_phi
        geo_B                 = self.geo_B
        D_km_interp           = interpolate.interp1d(geo_rho,geo_D)
        self.D_km_vals        = D_km_interp(rho_km_vals)
        rho_dot_interp        = interpolate.interp1d(geo_rho,geo_drho)
        self.rho_dot_kms_vals = rho_dot_interp(rho_km_vals)
        phi_deg_interp        = interpolate.interp1d(geo_rho,geo_phi)
        phi_deg_vals          = phi_deg_interp(rho_km_vals)
        B_deg_interp          = interpolate.interp1d(geo_rho,geo_B)
        B_deg_vals            = B_deg_interp(rho_km_vals)
        self.phi_rad_vals     = np.deg2rad(phi_deg_vals)
        self.B_rad_vals       = np.deg2rad(B_deg_vals)
        del geo_rho,geo_drho,geo_D,geo_phi,geo_B,D_km_interp,rho_dot_interp
        del phi_deg_vals,phi_deg_interp,B_deg_interp,B_deg_vals
