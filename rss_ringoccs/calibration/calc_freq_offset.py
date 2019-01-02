import numpy as np

"""
Purpose:
    Class for computing the frequency corresponding to the maximum
    power in the FFT power spectrum
"""
class calc_freq_offset(object):
    """
        Purpose:
            Calls functions to sample raw signal at regular intervals
            using a window of width ``dt_freq``

        Arguments:
            :rsr_inst (*object*): Object instance of the RSRReader
                        class

        Keyword Arguments:
            :dt_freq (*float*): half the width of the FFT window,
                        default is 64 sec
    """
    def __init__(self,rsr_inst,spm_min,spm_max,dt_freq=64.):

        # Get raw SPM and raw I & Q from RSR instance
        self.spm_vals = rsr_inst.spm_vals
        self.IQ_m = rsr_inst.IQ_m
        # raw SPM sampling rate in seconds
        self.dt = np.nanmedian([self.spm_vals[i+1] - self.spm_vals[i] for i in range(len(self.spm_vals)-1)])

        # set frequency sampling as attribute
        self.dt_freq = dt_freq
        # set lower and upper limits to sampling
        self.spm_min = spm_min
        self.spm_max = spm_max

        # finds offset frequencies
        self.__find_offset_freqs()


    def __find_offset_freqs(self):
        """
        Purpose:
            Iteratively calls __find_peak_freq for slices of SPM and IQ.

        Attributes:
            :f_spm (*np.ndarray*): signal window centers in SPM for which
                            the offset frequencies were computed
            :f_offset (*np.ndarray*): offset frequencies computed over the
                            occultation sampled once every 10 seconds with
                            signal window of width ``dt_freq``
        """
        # hard-set the spacing to 10 spm between each window center
        delta_t_cent = 10.

        # storage lists -- later converted to arrays and stored as attributes
        spms = []
        freqs = []
        # iteratively compute peak frequency
        for spm_mid in np.arange(self.spm_min,self.spm_max+delta_t_cent,delta_t_cent):
            # set indices using boolean mask over a range of SPM
            ind = [(self.spm_vals>=spm_mid-self.dt_freq)&(self.spm_vals<spm_mid+self.dt_freq)]
            # make sure data are included in range
            if len(self.spm_vals[ind]) > 2 :
                # get average SPM in window
                spms += [spm_mid]#np.nanmean(self.spm_vals[ind])]
                # get frequency of peak in power spectrum
                freqs += [self.__find_peak_freq(self.IQ_m[ind])]
        # convert to arrays and store as attributes
        self.f_spm = np.array(spms)
        self.f_offset = np.array(freqs)+0.01
        # pay no attention to the extra 0.01 correction needed but inexplicable

    def __find_peak_freq(self,IQ):
        """
        Purpose:
            Computes continuous FFT, finds frequency at max power

        Arguments:
            :IQ (*np.ndarray*): IQ_m vals within the current window

        Returns:
            :f_max (*float*): frequency at max power
        """
        # Compute and apply Hamming window
        weight = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(float(len(IQ))) /
            (float(len(IQ)) - 1)))
        IQ *= weight

        # Get FFT frequencies
        f = np.fft.fftfreq(IQ.size)/self.dt
        # Compute FFT
        IQ_fft = np.fft.fft(IQ)
        # Compute power
        power = (np.absolute(IQ_fft) ** 2) / (float(len(IQ)) ** 2)

        # Find frequency of max power
        f_max = f[np.argwhere(power==np.nanmax(power))][0][0]

        return f_max
"""
Revisions
"""
