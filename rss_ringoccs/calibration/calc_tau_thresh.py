"""

:Purpose:
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

:Dependencies:
    #. numpy
    #. matplotlib
    #. scipy

"""

import numpy as np
from scipy.signal import spectrogram

class calc_tau_thresh(object):
    """
    :Purpose:
        Compute threshold optical depth following

    Arguments
        :rsr_inst (*object*): object instance of the RSRReader class
        :geo_inst (*object*): object instance of the Geometry class
        :cal_inst (*object*): object instance of the Calibration class

    Keyword Arguments
        :res_km (*float*): Reconstruction resolution in km
        :Calpha (*float*): Constant for scaling Bandwidth/SNR ratio.
                        Default is 2.41 for 70% confidence
                        (see [MTR1986]_)

    Attributes
        :snr (*np.ndarray*): Signal-to-noise ratio SNR0 over the
                        entire occultation. This changes over the
                        occultation because the signal power
                        fluctuates.
        :tau_thresh (*np.ndarray*): threshold optical depth computed
                        using [MTR1986]_
        :spm_vals (*np.ndarray*): Observed event time array from cal_inst
        :rho_vals (*np.ndarray*): Ring intercept point array interpolated to
                                  spm_vals

    """
    def __init__(self,rsr_inst,geo_inst,cal_inst,
                res_km=1.0,Calpha=2.41):

        # find time-series sampling frequency and rate
        df = rsr_inst.sample_rate_khz * 1e3

        # resample rho and rho dot
        rho_km = np.interp(cal_inst.t_oet_spm_vals,geo_inst.t_oet_spm_vals,
                geo_inst.rho_km_vals)
        rho_dot_kms = np.interp(cal_inst.t_oet_spm_vals,geo_inst.t_oet_spm_vals,
                geo_inst.rho_dot_kms_vals)
        B_rad = np.deg2rad(np.interp(cal_inst.t_oet_spm_vals,
            geo_inst.t_oet_spm_vals,geo_inst.B_deg_vals))

        # set attributes
        self.spm_vals = cal_inst.t_oet_spm_vals
        self.rho_vals = rho_km

        # find noise power
        noise = self.find_noise(rsr_inst.spm_vals,rsr_inst.IQ_m,df)

        # find signal power
        signal = cal_inst.p_free_vals

        bandwidth = abs( rho_dot_kms / res_km )
        snr = signal/noise

        tau = -np.sin(abs(B_rad))*np.log(0.5*Calpha*bandwidth/snr)

        # compute SNR and set attribute
        self.snr = snr
        self.bandwidth = bandwidth
        self.mu = -np.sin(abs(B_rad))

        # compute threshold optical depth
        self.tau_thresh = tau

    def find_noise(self,spm,IQ,df):
        """
        Locate the additive receiver noise within the data set.
        This is done by computing a spectrogram of the raw
        complex signal, filtering out the spacecraft signal, and
        averaging over the frequency and time domains.

        Arguments
            :spm (*np.ndarray*): raw SPM in seconds
            :IQ (*np.ndarray*): measured complex signal
            :df (*float*): sampling frequency in Hz of the IQ

        Returns
            :p_noise (*np.ndarray*): noise power
        """
        # Spectrogram to filter out spacecraft signal
        n = int(1.024 * df) # number of elements based on frequency sampling
        freq,time,Sxx = spectrogram(IQ,df,nperseg=n,return_onesided=False)
        time += spm[0]
        # frequency filtering to include only the thermal receiver power,
        #   averaged over time
        Sxx_freq_filt = np.nanmean(Sxx[((freq>-450)&(freq<-200))|
            ((freq>200)&(freq<450)),:],0)
        # time filtering to include only the signal outside the occultation
        # by looking at just the first and last 1,000 seconds
        tmin = spm[0]+1e3
        tmax = spm[-1]-1e3
        p_noise = np.nanmean(Sxx_freq_filt[(time<tmin)|(time>tmax)])

        return p_noise
"""
Revision history
    Nov 29 2018 - sflury -- original
    Dec 14 2018 - sflury -- spectrogram method with filtering of space-
                            craft signal from FFT to get noise values
    Jan 11 2018 - sflury -- updated doc strings and method for
                            selecting the first and last thousand secs
                            of the occultation data set
    Nov 13 2020 - rmaguire -- Removed unused imports:
        import matplotlib.pyplot as plt
        from scipy.interpolate import splrep,splev
"""
