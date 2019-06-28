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
from scipy.signal import spectrogram

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
    def __init__(self,spm_vals,IQ,rho_km,rho_dot_kms,B_rad,p_params,
                ring_cent,ring_width,res_km=1.0,Calpha=2.49):

        # find time-series sampling frequency and rate
        df = round(float(len(spm_vals))/(spm_vals[-1]-spm_vals[0]))

        # set attributes
        self.spm_vals = spm_vals
        self.rho_vals = rho_km

        # find noise power
        noise = self.find_noise(spm_vals,IQ,df)

        # find signal power
        signal = np.polyval(p_params,self.rho_vals)

        bandwidth = abs( rho_dot_kms / res_km )
        snr = signal/noise

        # compute raw threshold optical depth
        tau = -np.sin(abs(B_rad))*np.log(0.5*Calpha*bandwidth/snr)
        # replace ring-region values (contaminated) with a fit from
        # surrounding computed threshold optical depth
        rmin = ring_cent-2*ring_width/np.mean(rho_dot_kms)
        rmax = ring_cent+2*ring_width/np.mean(rho_dot_kms)
        mask = [(spm_vals<rmin)|(spm_vals>rmax)]
        p = np.polyfit(spm_vals[mask],tau[mask],9)
        for i in range(len(spm_vals)):
            if spm_vals[i] >= rmin and spm_vals[i] <= rmax:
                tau[i] = np.polyval(p,spm_vals[i])

        # set snr attribute
        self.snr = snr
        # set noise attribute
        self.noise = noise
        # set bandwidth attribute
        self.bandwidth = bandwidth
        # compute threshold optical depth
        self.tau_thresh = tau

    def find_noise(self,spm,IQ,df):
        """
        Purpose:
            Locate the additive receiver noise within the data set.
            This is done by computing a spectrogram of the raw
            complex signal, filtering out the spacecraft signal, and
            averaging over the frequency and time domains.

        Arguments:
            :spm (*np.ndarray*): raw SPM in seconds
            :IQ (*np.ndarray*): measured complex signal
            :df (*float*): sampling frequency in Hz of the IQ

        Returns:
            :p_noise (*np.ndarray*): noise power
        """
        # Spectrogram to filter out spacecraft signal
        n = int(1.024 * df) # number of elements based on frequency sampling
        freq,time,Sxx = spectrogram(IQ,df,nperseg=n,return_onesided=False)
        time += spm[0]
        # frequency filtering to include only the thermal receiver power,
        #   averaged over time
        Sxx_freq_filt = np.nanmean(Sxx[((freq>-700)&(freq<-200))|
            ((freq>200)&(freq<700)),:],0)
        # average over entire occultation
        p_noise = np.nanmean(Sxx_freq_filt[(time<time[0]+10)|(time>time[-1]-10)])

        return p_noise
"""
Revision history
    Nov 29 2018 - sflury -- original
    Dec 14 2018 - sflury -- spectrogram method with filtering of space-
                            craft signal from FFT to get noise values
    Jan 11 2018 - sflury -- updated doc strings and method for
                            selecting the first and last thousand secs
                            of the occultation data set
"""
