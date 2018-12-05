"""
Purpose:
    Compute :math:`\\tau_{thresh}` for use as a proxy for maximum reliable value
    of optical depth within the diffraction-limited profile.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep,splev
from scipy.signal import spectrogram
import pdb

class calc_tau_thresh(object):
    """
    Purpose:
        Compute threshold optical depth following 

    Arguments:
        :rsr_inst (*object*): object instance of the RSRReader class
        :geo_inst (*object*): object instance of the Geometry class
        :freespace_spm (*np.ndarray*): locations of freespace where intrinsic
                                spacecraft signal is observed
        :IQ_c (*np.ndarray*): phase-corrected complex spacecraft signal
        :pnorm_fit (*np.ndarray*): polynomial fit to the freespace power,
                                resampled to raw SPM

    Keyword Arguments:
        :res_km (*float*): maximum threshold for noise power relative
                                to spacecraft signal power. Default is 10%.
        :Calpha (*float*): constant for scaling Bandwidth/SNR ratio. Default is
                                2.41 for 70% confidence

    Attributes:
        :pn_med (*float*): Median receiver "noise" power
        :pn_rms (*float*): RMS receiver "noise" power as proxy for uncertainty in
                                noise.
        :snr (*np.ndarray*): Signal-to-Noise Ratio over the entire occultation.
                                This changes over the occultation because the signal
                                power fluctuates.
        :tau_thresh (*np.ndarray*): threshold optical depth computed by

        .. math::

            \\tau_{thresh} = -\\sin(B)\\ln\\left(\\frac{C_{\\alpha} \\dot{\\rho}}{\\text{SNR}\\Delta R}\\right)

    """
    def __init__(self,rsr_inst,geo_inst,freespace_spm,IQ_c,pnorm_fit,
                res_km=1.0,Calpha=1.205,constant=False,raw=True):

        # get atmosphere location in observation
        atmo = geo_inst.atmos_occ_spm_vals

        ###
        ### ~~ RESAMPLE TO RAW SPM ~~
        ###

        # Compute spline coefficients relating SPM to rho and find rho
        rho_geo_spl_coef = splrep(geo_inst.t_oet_spm_vals, geo_inst.rho_km_vals)
        rho_vals = splev(rsr_inst.spm_vals, rho_geo_spl_coef)

        self.spm_vals = rsr_inst.spm_vals
        self.rho_vals = rho_vals

        # get bandwidth from spacecraft velocity relative to occultation event radius
        # and the data resolution, first resampling rho dot to raw
        rho_dot_geo_spl_coef = splrep(geo_inst.t_oet_spm_vals, geo_inst.rho_dot_kms_vals)
        rho_dot_vals = splev(rsr_inst.spm_vals,rho_dot_geo_spl_coef)
        bandwidth = abs(rho_dot_vals/res_km)

        # resample and convert ring elevation to radiians
        B_deg_geo_spl_coef = splrep(geo_inst.t_oet_spm_vals, geo_inst.B_deg_vals)
        B_deg_vals = splev(rsr_inst.spm_vals,B_deg_geo_spl_coef)
        B_rad = np.deg2rad(B_deg_vals)

        #self.spm_vals = rsr_inst.spm_vals
        #self.rho_vals = rho_vals

        ###
        ### ~~ FIND SIGNAL AND NOISE POWER ~~
        ###
        # noise
        noise = self.find_noise(rsr_inst.spm_vals,rsr_inst.IQ_m,pnorm_fit,freespace_spm)

        # signal
        # If constant SNR desired, use power from freespace regions
        if constant:
            signal = self.find_signal(spm_spec,spec,freespace_spm)
        # If SNR varies, use fit to freespace signal
        else:
            signal = pnorm_fit

        bandwidth = abs( rho_dot_vals / res_km )
        snr_raw = signal / noise
        tau_raw = np.sin(B_rad)*np.log(Calpha*bandwidth/snr_raw)
        self.snr = snr_raw
        self.tau_thresh = tau_raw

    def find_noise(self,spm,IQ,pnorm,freespace_spm):

        # boundaries of occultation
        tmin = np.nanmin(freespace_spm)
        tmax = np.nanmax(freespace_spm)

        # set maximum threshold for noise power relative to spacecraft signal
        p_sc_min = np.nanmin(pnorm[(spm>tmin)&(spm<tmax)])
        pnoise_thresh = 0.1 * p_sc_min
        power = abs(IQ**2)
        #
        noise = []
        for i in range(len(spm)):#spm,P in zip(spm_vals,power):
            if spm[i] < tmin or spm[i] > tmax :
                if i > 1 and i < len(spm)-2 :
                    if power[i] < pnoise_thresh and power[i-1] < pnoise_thresh:
                        noise += [power[i]]#[IQ[i]]#

        return np.nanmedian(noise)

    def find_noise_spec(self,f_spec,spm_spec,spec,pnorm_spec,freespace_spm):

        # boundaries of occultation
        tmin = np.nanmin(freespace_spm)
        tmax = np.nanmax(freespace_spm)
        pmin = np.nanmin(pnorm_spec[(spm_spec>tmin)&(spm_spec<tmax)])

        # Where noise is located in each spectrum
        noise_ind_hz = (((f_spec >= -400) &
            (f_spec <= -200)) |
            ((f_spec >= 200) &
             (f_spec <= 400)))

        # Noise and spectrum signal for each spectrum
        noise_vals = np.mean(spec[noise_ind_hz, :], 0)

        # Indices to use for thermal noise before occultation begins
        # from the gap finding routine
        noise_ind_spm = np.argwhere( (spm_spec <= tmin) | (spm_spec >= tmax) )
        noise_ind = np.reshape(np.concatenate(noise_ind_spm), -1)

        # clip based on power
        pclip = np.argwhere(noise_vals[noise_ind]<0.1*pmin)
        return np.nanmedian(noise_vals[noise_ind][pclip])



    def find_signal(self,spm_spec,spec,freespace_spm):

        # Indices to estimate signal using freespace regions
        signal_ind = []
        for i in range(len(freespace_spm)):
            signal_ind.append(np.argwhere(
                (spm_spec >= freespace_spm[i][0]) &
                (spm_spec <= freespace_spm[i][1])))
        signal_ind = np.reshape(np.concatenate(signal_ind), -1)

        return np.nanmedian(spec[signal_ind])

"""
Revision history
    Nov 29 2018 - sflury -- original
"""
