"""
    Module Name:
        diffcorr
    Purpose:
        Provide functions and classes that aid in the process of
        diffraction Correction / Fresnel Inversion. Additional
        functions for the purpose of forward modelling of
        reconstructed data and diffraction modelling exists.
        Several low-level functions that perform error checks for
        the main functions are also included.
    
    Window (Taper) Functions:
        rect................Rectangular window.
        coss................Squared cossine window.
        kb25................Kaiser-Bessel 2.5 window.
        kb35................Kaiser-Bessel 3.5 window.
        kbmd................Modified Kaiser-Bessel 2.5 window.

    Error Checks:
        check_boole.........Checks for Boolean input.
        check_ternary.......Checks for Ternary input.
        check_pos_real......Checks if an input is positive and real.
        check_real..........Checks if input is real (Number/Array).
        check_complex.......Checks if input is complex.
    
    Special Functions:
        fresnel_sin.........The Fresnel sine integral.
        fresnel_cos.........The Fresnel cosine integral.
        sq_well_solve.......Diffraction pattern through square well.

    Mathematical Functions:
        compute_norm_eq.....Computes the normalized equivalent width.
        get_norm_eq.........Quickly retrieve pre-computed normalized
                            equivalent widths from strings with the
                            name of common window functions.
        resolution_inverse..Computes the inverse of the function
                            y = x/(exp(-x)+x-1)
        power_func..........Compute power from complex transmittance.
        phase_func..........Compute phase from complex transmittance.
        tau_func............Compute normalized optical depth from the
                            complex transmittance.
        wker................Computes a weighted kernel function.
        freq_wav............Convert frequency to wavelength, and
                            vice-versa. Kilometers or Herts only.
        fresnel_scale.......Compute the Fresnel scale.
    
    Miscellaneous Functions:
        get_range_request...Computes an array of the form [a,b] from
                            a given array, list, or from a set of
                            allowed strings.
        get_range_actual....Given an array of numbers (usually the
                            radial range of the data), a range
                            request, and a window width, compute the
                            allowed range of processing.
"""

# Import dependencies for the diffcorrpy module
import spiceypy as spice
import numpy as np
import pandas as pd
import time
import sys
import subprocess
from .misc_funcs import *
from scipy.special import lambertw
from scipy.special import iv
from scipy.special import erf
from scipy import interpolate

class rec_data(object):
    """
        Class: rec_data
        Purpose: This object contains all of input and output variables from
        diffraction reconstruction, including geometry, data, and reconstructed
        data.
        Attributes:
            rho_km_vals (np.ndarray): Ring intercept radius values, in km,
                at final spacing specified in dr_km
            p_norm_vals (np.ndarray): Power normalized to 1. This is the diffraction
                pattern that is input to the Fresnel inversion step
            phase_rad_vals (np.ndarray): Phase of the complex signal, in radians.
                This is the other part of the diffraction pattern that is input to
                the Fresnel Inversion step
            B_rad_vals (np.ndarray): Ring opening angle of Saturn
            D_km_vals (np.ndarray): Spacecraft-RIP (Ring Intercept Point) distance
            f_sky_hz_vals (np.ndarray): Predicted sky frequency, in Hz.
            phi_rad_vals (np.ndarray): Observed ring azimuth, in radians
            rho_dot_kms_vals (np.ndarray): Ring intercept point velocity
	"""

    def __init__(self, NormDiff,res,wtype,bfac=True):
        """
            Arguments:
                NormDiff: Instance of the class NormDiff containing 
        """

        self.res                = None
        self.wtype              = None
        self.rho_km_vals        = None
        self.p_norm_vals        = None
        self.phase_rad_vals     = None
        self.b_rad_vals         = None
        self.d_km_vals          = None
        self.f_sky_hz_vals      = None
        self.phi_rad_vals       = None
        self.rho_dot_kms_vals   = None
        self.T_hat_vals         = None
        self.F_km_vals          = None
        self.w_km_vals          = None
        self.mu_vals            = None
        self.lambda_sky_km_vals = None
        self.dx_km              = None
        self.norm_eq            = None
        self.res                = res
        self.wtype              = wtype

        self.__set_attributes(NormDiff)
        self.__get_occ_type()

        n_rho = np.size(self.rho_km_vals)
        self.__error_check(n_rho)

        self.lambda_sky_km_vals = freq_wav(self.f_sky_hz_vals)
        self.dx_km              = self.rho_km_vals[1] - self.rho_km_vals[0]
        self.T_hat_vals         = np.sqrt(np.abs(self.p_norm_vals)) * \
            np.exp(-1j * self.phase_rad_vals)
        self.F_km_vals          = fresnel_scale(self.lambda_sky_km_vals,
            self.d_km_vals,self.phi_rad_vals,self.b_rad_vals)
        self.norm_eq            = get_norm_eq(self.wtype)
        self.w_km_vals          = window_width(self.res,self.norm_eq,
            self.f_sky_hz_vals,self.F_km_vals,self.rho_dot_kms_vals,bfac=bfac)
        self.mu_vals            = np.sin(np.abs(self.b_rad_vals))
    
    def __set_attributes(self, NormDiff):
        self.rho_km_vals      = NormDiff.rho_km_vals
        self.p_norm_vals      = NormDiff.p_norm_vals
        self.phase_rad_vals   = NormDiff.phase_rad_vals
        self.b_rad_vals		  = NormDiff.b_rad_vals
        self.d_km_vals		  = NormDiff.d_km_vals
        self.f_sky_hz_vals    = NormDiff.f_sky_hz_vals
        self.phi_rad_vals     = NormDiff.phi_rad_vals
        self.rho_dot_kms_vals = NormDiff.rho_dot_kms_vals
    
    def __error_check(self,N_RHO):
        error_code = 0
        if (np.size(self.phase_rad_vals) != N_RHO):
            print("Bad Input: len(phase_rad_vals) != len(rho_km_vals)")
            error_code += 1
        if (np.size(self.p_norm_vals) != N_RHO):
            print("Bad Input: len(p_norm_vals) != len(rho_km_vals)")
            error_code += 2
        if (np.size(self.phi_rad_vals) != N_RHO):
            print("Bad Input: len(phi_rad_vals) != len(rho_km_vals)")
            error_code += 4
        if (np.size(self.b_rad_vals) != N_RHO):
            print("Bad Input: len(B_rad_vals) != (rho_km_vals)")
            error_code += 8
        if (np.size(self.f_sky_hz_vals) != N_RHO):
            print("Bad Input: len(f_sky_Hz_vals) != len(rho_km_vals)")
            error_code += 16
        if (np.size(self.d_km_vals) != N_RHO):
            print("Bad Input: len(D_km_vals) != len(rho_km_vals)")
            error_code += 32
        if (np.size(self.rho_dot_kms_vals) != N_RHO):
            print("Bad Input: len(D_km_vals) != len(rho_km_vals)")
            error_code += 64
        if (error_code != 0):
            sys.exit("Bad input. Read error message.")

    def __get_occ_type(self):
        drho = [np.min(self.rho_dot_kms_vals),np.max(self.rho_dot_kms_vals)]
        dx   = self.rho_km_vals[1]-self.rho_km_vals[0]
        if (drho[0] < 0) and (drho[1] > 0):
            sys.exit("drho/dt had positive and negative values.")
        if (dx > 0) and (drho[1] < 0):
            self.rho_dot_kms_vals=np.abs(self.rho_dot_kms_vals)
        if (dx < 0):
            self.rho_km_vals      = self.rho_km_vals[::-1]
            self.phase_rad_vals   = self.phase_rad_vals[::-1]
            self.p_norm_vals      = self.p_norm_vals[::-1]
            self.phi_rad_vals     = self.phi_rad_vals[::-1]
            self.b_rad_vals       = self.b_rad_vals[::-1]
            self.f_sky_hz_vals    = self.f_sky_hz_vals[::-1]
            self.d_km_vals        = self.d_km_vals[::-1]
            self.rho_dot_kms_vals = np.abs(self.rho_dot_kms_vals[::-1])

class diffraction_correction(object):

    def __init__(self,dat,res,rng="all",wtype="kb25",dir="",rev="",
        fwd=False,norm=True,fast=True,verbose=True):
        t1       = time.time()

        self.res                = None
        self.wtype              = None
        self.rng                = None
        self.rho_km_vals        = None
        self.p_norm_vals        = None
        self.phase_rad_vals     = None
        self.b_rad_vals         = None
        self.d_km_vals          = None
        self.f_sky_hz_vals      = None
        self.phi_rad_vals       = None
        self.rho_dot_kms_vals   = None
        self.T_hat_vals         = None
        self.F_km_vals          = None
        self.w_km_vals          = None
        self.mu_vals            = None
        self.lambda_sky_km_vals = None
        self.dx_km              = None
        self.norm_eq            = None
        self.n_used             = None
        self.start              = None
        self.T_vals             = None
        self.power_vals         = None
        self.tau_vals           = None
        self.phase_vals         = None
        self.p_norm_fwd_vals    = None
        self.T_hat_fwd_vals     = None
        self.phase_fwd_vals     = None
        self.norm               = None
        self.fast               = None
        self.fwd                = None

        self.res   = res
        self.wtype = wtype
        self.rng   = rng
        self.norm  = norm
        self.fast  = fast
        self.fwd   = fwd

        recdata    = rec_data(dat, res, wtype)

        self.__set_dc_attributes(recdata)
        self.__compute_dc_attributes(rng)

        self.T_vals = self.__fresnel_inversion_fast()
        
        if fwd:
            self.T_hat_fwd_vals=fresnel_forward(self.rho_km_vals,\
            self.F_km_vals,self.phi_rad_vals,self.b_rad_vals,\
            self.d_km_vals,self.T_vals,self.lambda_sky_km_vals,\
            self.w_km_vals,self.dx_km,self.wtype,self.start,self.n_used,\
            Normalize=True)
        self.power_vals = power_func(self.T_vals)
        self.phase_vals = phase_func(self.T_vals)
        self.tau_vals   = tau_func(self.T_vals,self.mu_vals)

        self.__trim_attributes(fwd)

        t2 = time.time()
        print("Computation Time: ",t2-t1)

    def __rect(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        w_func = np.zeros(nw_pts) + 1.0
        return w_func

    def __coss(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
        w_func = np.cos(np.pi * x / w_in)**2
        return w_func

    def __kb25(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.5*np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kb35(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 3.5 * np.pi
        w_func = iv(0.0,alpha * np.sqrt((1.0 - (2.0 * x / w_in)**2)))/iv(0.0,alpha)
        return w_func

    def __kbmd(w_in, dx):
        nw_pts = int(2 * np.floor(w_in / (2.0 * dx)) + 1)
        x      = (np.array((range(nw_pts))) - ((nw_pts - 1) / 2.0)) * dx
        alpha  = 2.5 * np.pi
        w_func = (iv(0.0,alpha*np.sqrt((1.0-(2.0*x/w_in)**2)))-1)/(iv(0.0,alpha)-1)
        return w_func
    
    __func_dict = {
        "__rect" : {"func" : __rect, "normeq" : 1.00000000},
        "__coss" : {"func" : __coss, "normeq" : 1.50000000},
        "__kb25" : {"func" : __kb25, "normeq" : 1.65191895},
        "__kb35" : {"func" : __kb35, "normeq" : 1.92844639},
        "__kbmd" : {"func" : __kbmd, "normeq" : 1.65994218}
        }

    def __set_dc_attributes(self,recdata):
        self.res                = recdata.res
        self.wtype              = recdata.wtype
        self.rho_km_vals        = recdata.rho_km_vals
        self.p_norm_vals        = recdata.p_norm_vals
        self.phase_rad_vals     = recdata.phase_rad_vals
        self.b_rad_vals         = recdata.b_rad_vals
        self.d_km_vals          = recdata.d_km_vals
        self.f_sky_hz_vals      = recdata.f_sky_hz_vals
        self.phi_rad_vals       = recdata.phi_rad_vals
        self.rho_dot_kms_vals   = recdata.rho_dot_kms_vals
        self.T_hat_vals         = recdata.T_hat_vals
        self.F_km_vals          = recdata.F_km_vals
        self.w_km_vals          = recdata.w_km_vals
        self.mu_vals            = recdata.mu_vals
        self.lambda_sky_km_vals = recdata.lambda_sky_km_vals
        self.dx_km              = recdata.dx_km
        self.norm_eq            = recdata.norm_eq
    
    def __compute_dc_attributes(self,rng):
        self.rng               = get_range_request(rng)
        self.start,self.n_used = get_range_actual(self.rho_km_vals,self.rng,
            self.w_km_vals)
    
    def __trim_attributes(self,fwd):
        start  = self.start
        n_used = self.n_used
        crange = np.arange(n_used)+start

        self.rho_km_vals         = self.rho_km_vals[crange]
        self.p_norm_vals         = self.p_norm_vals[crange]
        self.phase_rad_vals      = self.phase_rad_vals[crange]
        self.b_rad_vals          = self.b_rad_vals[crange]
        self.d_km_vals           = self.d_km_vals[crange]
        self.f_sky_hz_vals       = self.f_sky_hz_vals[crange]
        self.phi_rad_vals        = self.phi_rad_vals[crange]
        self.rho_dot_kms_vals    = self.rho_dot_kms_vals[crange]
        self.T_hat_vals          = self.T_hat_vals[crange]
        self.F_km_vals           = self.F_km_vals[crange]
        self.w_km_vals           = self.w_km_vals[crange]
        self.mu_vals             = self.mu_vals[crange]
        self.lambda_sky_km_vals  = self.lambda_sky_km_vals[crange]
        self.T_vals              = self.T_vals[crange]
        self.power_vals          = self.power_vals[crange]
        self.tau_vals            = self.tau_vals[crange]
        self.phase_vals          = self.phase_vals[crange]
        if fwd:
            self.p_norm_fwd_vals = self.p_norm_fwd_vals[crange]
            self.T_hat_fwd_vals  = self.T_hat_fwd_vals[crange]
            self.phase_fwd_vals  = self.phase_fwd_vals[crange]
    
    def __trim_inputs(self):
        start  = self.start
        n_used = self.n_used
        w   = np.max(self.w_km_vals[start:start+n_used])
        nst = np.min((self.rho_km_vals>=(self.rho_km_vals[start]-w)).nonzero())
        nen = np.max((self.rho_km_vals<=(self.rho_km_vals[start]+w)).nonzero())
        nreq   = nen - nst
        crange = np.array(range(0,nst))+nreq+1
        print(np.size(self.rho_km_vals))
        self.rho_km_vals         = self.rho_km_vals[crange]
        print(np.size(self.rho_km_vals))
        self.start               = start-nreq
        self.p_norm_vals         = self.p_norm_vals[crange]
        self.phase_rad_vals      = self.phase_rad_vals[crange]
        self.b_rad_vals          = self.b_rad_vals[crange]
        self.d_km_vals           = self.d_km_vals[crange]
        self.f_sky_hz_vals       = self.f_sky_hz_vals[crange]
        self.phi_rad_vals        = self.phi_rad_vals[crange]
        self.rho_dot_kms_vals    = self.rho_dot_kms_vals[crange]
        self.T_hat_vals          = self.T_hat_vals[crange]
        self.F_km_vals           = self.F_km_vals[crange]
        self.w_km_vals           = self.w_km_vals[crange]
        self.mu_vals             = self.mu_vals[crange]
        self.lambda_sky_km_vals  = self.lambda_sky_km_vals[crange]

    def __fresinv(self,T_hat,ker,dx,f_scale):
        T = np.sum(ker * T_hat) * dx * (1.0+1.0j) / (2.0 * f_scale)
        return T

    def __psifacfast(self,r,r0,cb,cp0,sp0):
        factor  = ((cb*cb) * cp0 * sp0 / (1.0 - (cb*cb) * (sp0*sp0))) * (r - r0) / r0
        return factor

    def __psifast(self,r,r0,d,cb,cp,sp,cp0,sp0):
        xi   = (cb / d) * (r0*cp0 - r*cp)
        eta  = ((r0*r0) + (r*r) - 2.0 * r * r0 * (sp*sp0 + cp*cp0)) / (d*d)
        psi_vals   = np.sqrt(1.0 + 2.0 * xi + eta) - (1.0 + xi)
        return psi_vals

    def __normalize(self,r,w_func,f_scale):
        x         = r-np.mean(r)
        drho      = r[1]-r[0]
        psi       = (np.pi / 2.0) * ((x / f_scale)**2)
        ker       = np.exp(1j * psi)
        T1        = np.abs(np.sum(w_func * ker) * drho)
        norm_fact = np.sqrt(2.0) * f_scale / T1
        return norm_fact

    def __fresnel_inversion_fast(self):
        # Retrieve variables.
        w_vals       = self.w_km_vals
        rho_vals     = self.rho_km_vals
        phi_rad_vals = self.phi_rad_vals
        d_vals       = self.d_km_vals
        b_rad_vals   = self.b_rad_vals
        lambda_vals  = self.lambda_sky_km_vals
        T_hat_vals   = self.T_hat_vals
        F_vals       = self.F_km_vals
        norm         = self.norm
        wtype        = "%s%s" % ("__",self.wtype)
        start        = self.start
        n_used       = self.n_used
        dx           = self.dx_km
        # Compute necessary variables.
        kD_vals   = 2.0 * np.pi * d_vals / lambda_vals
        cosb      = np.cos(b_rad_vals)
        cosphi0   = np.cos(phi_rad_vals)
        sinphi0   = np.sin(phi_rad_vals)
        dsq       = d_vals*d_vals
        rsq       = rho_vals*rho_vals
        # Define functions
        fw        = self.__func_dict[wtype]["func"]
        psifac    = self.__psifacfast
        finv      = self.__fresinv
        psif      = self.__psifast
        nrm       = self.__normalize
        # Calculate the corrected complex amplitude, point by point
        T_vals    = T_hat_vals * 0.0
        w_init    = w_vals[start]
        w_func    = fw(w_init,dx)
        nw        = np.size(w_func)
        phi_s_rad1 = phi_rad_vals[start]
        for i in np.arange(n_used):
            center = start+i
            r0     = rho_vals[center]
            r02    = rsq[center]
            d      = d_vals[center]
            d2     = dsq[center]
            cb     = cosb[center]
            cp0    = cosphi0[center]
            sp0    = sinphi0[center]
            kD     = kD_vals[center]
            w      = w_vals[center]
            if (w_init - w>= 2.0*dx):
                w_init     = w
                w_func     = fw(w,dx)
                nw         = np.size(w_func)
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                dphi_s_rad = psifac(r,r0,cb,cp0,sp0)
                phi_s_rad  = phi_rad_vals[center] - dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
            else:
                crange     = np.arange(int(center-(nw-1)/2),int(1+center+(nw-1)/2))
                r          = rho_vals[crange]
                r2         = rsq[crange]
                phi_s_rad  = phi_s_rad1
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)

            phi_s_rad1 = phi_s_rad
            loop = 0

            # Perform Newton-Raphson on phi.
            while (np.max(np.abs(dphi_s_rad)) > 1.e-8):
                xi         = (cb / d) * (r0*cp0 - r*cp)
                eta        = (r02 + r2 - 2.0*r*r0*(sp*sp0 + cp*cp0)) / d2
                v1         = r * cb * sp / d
                v2         = 2.0 * r * r0 * (sp*cp0 - sp0*cp) / (d2)
                v3         = 2.0*v1 + v2
                v4         = np.sqrt(1.0 + 2.0*xi + eta)
                v5         = r * cb * cp / d
                v6         = 2.0 * r * r0 * (sp*sp0 + cp*cp0) / d2
                dphia      = (2.0*v5 + v6)/(2.0 * v4)
                dphib      = v5 + (2.0*v1 + v2)*(2.0*v1 + v2)/(4.0*(v4*v4*v4))
                psi_d1     = v3 / (2.0 * v4) - v1
                psi_d2     = dphia - dphib
                dphi_s_rad = -psi_d1 / psi_d2
                phi_s_rad += dphi_s_rad
                cp         = np.cos(phi_s_rad)
                sp         = np.sin(phi_s_rad)
                loop      += 1
                if loop > 5:
                    break
            
            # Compute psi and then compute the forward model.
            psi_vals = kD * psif(r,r0,d,cb,cp,sp,cp0,sp0)
            ker      = w_func*np.exp(-1j*psi_vals)
            T_hat    = T_hat_vals[crange]
            F        = F_vals[center]
            T_vals[center] = finv(T_hat,ker,dx,F)
            if norm:
                T_vals[center] *= nrm(r,w_func,F)
            #print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
            #% (i,n_used,nw,loop),end="\r")
        #print("Pt: %d  Tot: %d  Width: %d  Psi Iters: %d  Fast Inversion" \
        #% (i,n_used,nw,loop))
        return T_vals

class compare_tau(object):
    def __init__(self,geodata,caldata,dlpdata,taudata,occ,res,rng='all'):
        data             = extract_csv_data(geodata,caldata,dlpdata,occ)
        tau_dat          = self.__read_tau_file(taudata)
        rec              = diffraction_correction(data,res,rng=rng)
        rho_km_vals      = rec.rho_km_vals
        tr               = tau_dat[...,0]
        tt               = tau_dat[...,5]
        tbdeg            = tau_dat[...,11]
        rmin             = np.min(tr)
        rmax             = np.max(tr)
        rstart           = int(np.min((rho_km_vals-rmin>=0).nonzero()))
        rfin             = int(np.max((rmax-rho_km_vals>=0).nonzero()))
        rho_km_vals      = rho_km_vals[rstart:rfin]
        tbrad            = tbdeg*spice.rpd()
        tm               = np.sin(np.abs(tbrad))
        tt_interp        = interpolate.interp1d(tr,tt,kind='linear')
        tau_tau          = tt_interp(rho_km_vals)
        tm_intperp       = interpolate.interp1d(tr,tm,kind='linear')
        tmu              = tm_intperp(rho_km_vals)
        tau_power        = np.exp(-tau_tau/tmu)
        self.rho_km_vals = rho_km_vals
        power_vals       = rec.power_vals
        self.power_vals  = power_vals[rstart:rfin]
        tau_vals         = rec.tau_vals
        self.tau_vals    = tau_vals[rstart:rfin]
        self.tau_power   = tau_power
        self.tau_tau     = tau_tau

    def __read_tau_file(self,taudata):
        dft = np.genfromtxt(taudata, delimiter=',')
        return dft