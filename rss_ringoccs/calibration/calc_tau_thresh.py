"""
Purpose:
    Compute :math:`\\tau_{thresh}` for use as a proxy for maximum reliable value
    of optical depth within the diffraction-limited or diffraction-reconstructed
    profile. This follows [MTR1986]_ Equations 24 and 26, which define,
    respectively, the power of the thermal noise

        .. math::
            \\hat{P}_N = \\frac{\\dot{\\rho}_0}{\\mathrm{SNR}_0 \\Delta\\rho_0}

    and the threshold optical depth

        .. math::
            \\tau_{thresh} = -\\sin(B)\\ln\\left(\\frac{1}{2}C_{\\alpha}\\hat{P}_N\\right)

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
        :pnorm_fit (*np.ndarray*): polynomial fit to the freespace power,
                                resampled to raw SPM

    Keyword Arguments:
        :res_km (*float*):
        :Calpha (*float*): constant for scaling Bandwidth/SNR ratio. Default is
                                2.41 for 70% confidence (see [MTR1986]_)

    Attributes:
        :snr (*np.ndarray*): Signal-to-noise ratio SNR0 over the entire occultation.
                                This changes over the occultation because the signal
                                power fluctuates.
        :tau_thresh (*np.ndarray*): threshold optical depth computed using [MTR1986]_

    """
    def __init__(self,rsr_inst,geo_inst,freespace_spm,pnorm_fit,
                res_km=1.0,Calpha=2.41,constant=False):

        # get atmosphere location in observation
        atmo = geo_inst.atmos_occ_spm_vals

        ###
        ### ~~ COMPUTE SPECTROGRAM ~~
        ###
        # Points per FFT for spectrogram
        fs = rsr_inst.sample_rate_khz*1000
        nperseg = int(fs*1.024)

        # computing the spectrogram
        f_spec, t_spec, spec = spectrogram(rsr_inst.IQ_m,fs,
            nperseg=int(fs*1.024), return_onesided=False)

        ###
        ### ~~ RESAMPLE TO SPECTROGRAM SPM ~~
        ###
        # SPM and various geometry/calibration parameters for spectra
        spm_spec = rsr_inst.spm_vals[0] + t_spec
        self.spm_vals = spm_spec
        rho_km_spec = np.interp(spm_spec, geo_inst.t_oet_spm_vals, geo_inst.rho_km_vals)
        self.rho_vals = rho_km_spec
        rho_dot_kms_spec = np.interp(spm_spec, geo_inst.t_oet_spm_vals, geo_inst.rho_dot_kms_vals)
        B_deg_spec = np.interp(spm_spec, geo_inst.t_oet_spm_vals, geo_inst.B_deg_vals)
        B_rad = np.deg2rad(B_deg_spec)
        pnorm_spec = np.interp(spm_spec, rsr_inst.spm_vals, pnorm_fit)

        ###
        ### ~~ FIND SIGNAL AND NOISE POWER ~~
        ###
        # noise
        noise = self.find_noise_spec(f_spec,spm_spec,spec,pnorm_spec,freespace_spm)
        # signal
        # If constant SNR desired, use power from freespace regions
        if constant:
            signal = self.find_signal(spm_spec,spec,freespace_spm)
        # If SNR varies, use fit to freespace signal
        else:
            signal = np.interp(spm_spec, rsr_inst.spm_vals, pnorm_fit)
        #
        bandwidth = abs( rho_dot_kms_spec / res_km )
        snr_spec = signal/noise
        tau_spec = np.sin(B_rad)*np.log(0.5*Calpha*bandwidth/snr_spec)

        # compute SNR and set attribute
        self.snr = snr_spec

        ###
        ### ~~ COMPUTE THRESHOLD OPTICAL DEPTH ~~
        ###
        self.tau_thresh = tau_spec

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
            (f_spec <= -100)) |
            ((f_spec >= 100) &
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
