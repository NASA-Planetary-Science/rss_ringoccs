"""
Purpose:
    Compute :math:`\\tau_{thresh}` for use as a proxy for maximum
    reliable value of optical depth within the diffraction-limited or
    diffraction-reconstructed profile. This follows [MTR1986]_
    Equations 24 and 26, which define, respectively, the power of the
    thermal noise

    .. math::
        \\hat{P}_N = \\frac{\\dot{\\rho}_0}
        {\\mathrm{SNR}_0 \\Delta\\rho_0}

    and the threshold optical depth

    .. math::
        \\tau_{thresh} = -\\sin(B)
        \\ln\\left(\\frac{1}{2}C_{\\alpha}\\hat{P}_N\\right)

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
        :freespace_spm (*np.ndarray*): locations of freespace where
                        intrinsic spacecraft signal is observed
        :pnorm_fit (*np.ndarray*): polynomial fit to the freespace
                        power, resampled to raw SPM

    Keyword Arguments:
        :res_km (*float*):
        :Calpha (*float*): constant for scaling Bandwidth/SNR ratio.
                        Default is 2.41 for 70% confidence
                        (see [MTR1986]_)

    Attributes:
        :snr (*np.ndarray*): Signal-to-noise ratio SNR0 over the
                        entire occultation. This changes over the
                        occultation because the signal power
                        fluctuates.
        :tau_thresh (*np.ndarray*): threshold optical depth computed
                        using [MTR1986]_

    """
    def __init__(self,rsr_inst,geo_inst,cal_inst,
                res_km=1.0,Calpha=2.41,constant=False):

        # find time-series sampling frequency and rate
        dt = np.nanmedian([rsr_inst.spm_vals[i+1]-rsr_inst.spm_vals[i] for i in range(len(rsr_inst.spm_vals)-1)])
        #cal_inst.t_oet_spm_vals[1]-cal_inst.t_oet_spm_vals[0]
        df = rsr_inst.sample_rate_khz * 1e3

        # resample rho and rho dot
        rho_km = np.interp(cal_inst.t_oet_spm_vals,geo_inst.t_oet_spm_vals,geo_inst.rho_km_vals)
        rho_dot_kms = np.interp(cal_inst.t_oet_spm_vals,geo_inst.t_oet_spm_vals,geo_inst.rho_dot_kms_vals)
        B_rad = np.deg2rad(np.interp(cal_inst.t_oet_spm_vals,geo_inst.t_oet_spm_vals,geo_inst.B_deg_vals))

        # set attributes
        self.spm_vals = cal_inst.t_oet_spm_vals
        self.rho_vals = rho_km

        ###
        ### ~~ FIND SIGNAL AND NOISE POWER ~~
        ###
        # noise
        #noise = self.find_noise_spec(f_spec,spm_spec,spec,pnorm,geo_inst.freespace_spm)
        noise = self.find_noise(rsr_inst.spm_vals,rsr_inst.IQ_m,dt,df)
        # signal
        #signal = np.interp(spm_spec, cal_inst.t_oet_spm_vals, cal_inst.p_free_vals)
        signal = cal_inst.p_free_vals
        #
        bandwidth = abs( rho_dot_kms / res_km )
        snr = signal/noise

        tau = -np.sin(abs(B_rad))*np.log(0.5*Calpha*bandwidth/snr)

        # compute SNR and set attribute
        self.snr = snr

        ###
        ### ~~ COMPUTE THRESHOLD OPTICAL DEPTH ~~
        ###
        self.tau_thresh = tau

    def find_noise(self,spm,IQ,dt,df):
        """
        Purpose:
            Locate the additive receiver noise within the data set.
            This is done by computing a spectrogram of the raw
            complex signal, filtering out the spacecraft signal, and
            averaging over the frequency and time domains.

        Arguments:
            :spm (*np.ndarray*): raw SPM in seconds
            :IQ (*np.ndarray*): measured complex signal
            :dt (*float*): sampling spacing in sec of the raw SPM for
                        use in determining the indices of the first
                        and last 1,000 seconds of the occultation in
                        order to determine the location of receiver
                        noise within the data set
            :df (*float*): sampling frequency in Hz of the IQ
        """
        # Spectrogram to filter out spacecraft signal
        n = int(1.024 * df) # number of elements based on frequency sampling
        freq,time,Sxx = spectrogram(IQ,df,nperseg=n,return_onesided=False)
        time += spm[0]
        # frequency filtering to include only the thermal receiver power, averaged over time
        Sxx_freq_filt = np.nanmean(Sxx[((freq>-450)&(freq<-200))|((freq>200)&(freq<450)),:],0)
        # time filtering to include only the signal outside the occultation
        # by looking at just the first and last 1,000 seconds
        tmin = spm[int(1e3/dt)]
        tmax = spm[-int(1e3/dt)]
        p_noise = np.nanmean(Sxx_freq_filt[(time<tmin)|(time>tmax)])

        return p_noise
"""
Revision history
    Nov 29 2018 - sflury -- original
"""
