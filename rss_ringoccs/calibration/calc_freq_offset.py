"""

:Purpose:
    Class for computing the frequency corresponding to the maximum
    power in the FFT power spectrum

:Dependencies:
    #. numpy
"""

import numpy as np

class calc_freq_offset(object):
    """
    :Purpose:
        Calls functions to sample raw signal at regular intervals
        using a window of width ``dt_freq``

    Arguments
        :rsr_inst (*object*): Object instance of the RSRReader class
        :spm_min (*float*): Minimum observed event time for sampling
        :spm_max (*float): Maximum observed event time for sampling

    Keyword Arguments
        :dt_freq (*float*): half the width of the FFT window

    Attributes
        :spm_vals (*np.ndarray*): Observed event time at full sampling
        :IQ_m (*np.ndarray*): Uncorrected real and imaginary components of signal
        :dt (*float*): Raw time sampling from spm_vals
        :dt_freq (*float*): half the width of the FFT window
        :spm_min (*float*): Minimum time for sampling
        :spm_max (*float*): Maximum time for sampling
        :f_spm (*np.ndarray*): Observed event time for frequency offset
        :f_offset (*np.ndarray*): Frequency offset, or frequency at max power

    """
    def __init__(self,rsr_inst,spm_min,spm_max,dt_freq=2.):

        # Get raw SPM and raw I & Q from RSR instance
        self.spm_vals = rsr_inst.spm_vals
        self.IQ_m = rsr_inst.IQ_m
        # raw SPM sampling rate in seconds
        self.dt = self.spm_vals[1]-self.spm_vals[0]

        # set frequency sampling as attribute
        self.dt_freq = dt_freq
        # set lower and upper limits to sampling
        self.spm_min = spm_min
        self.spm_max = spm_max

        # finds offset frequencies
        self.__find_offset_freqs()


    def __find_offset_freqs(self):
        """
        Iteratively calls __find_peak_freq for slices of SPM and IQ
        and sets f_spm and f_offset attributes.
        """
        # hard-set the spacing to 10 spm between each window center
        delta_t_cent = 10.

        # storage lists -- later converted to arrays and stored as attributes
        spms = []
        freqs = []
        # iteratively compute peak frequency
        for spm_mid in np.arange(self.spm_min,self.spm_max+delta_t_cent,
                delta_t_cent):
            # set indices using boolean mask over a range of SPM
            ind = [(self.spm_vals>=spm_mid-self.dt_freq)&
                    (self.spm_vals<spm_mid+self.dt_freq)]
            # make sure data are included in range
            if len(self.spm_vals[ind]) > 2 :
                # get average SPM in window
                spms += [spm_mid]
                # get frequency of peak in power spectrum
                freqs += [self.__find_peak_freq(self.spm_vals[ind],
                                                self.IQ_m[ind])]
        # convert to arrays and store as attributes
        self.f_spm = np.array(spms)
        self.f_offset = np.array(freqs)

    def __find_peak_freq(self,time,IQ,df=0.001,hwid=0.2):
        """
        Computes continuous FFT, finds frequency at max power

        Arguments
            :IQ (*np.ndarray*): IQ_m vals within the current window

        Returns
            :f_max (*float*): frequency at max power
        """
        # Compute and apply Hamming window
        weight = 0.5 * (1.0 - np.cos(2.0 * np.pi *
                    np.arange(float(len(IQ))) / (float(len(IQ)) - 1)))
        IQ *= weight

        # Get FFT frequencies
        f = np.fft.fftfreq(IQ.size,d=self.dt)
        # Compute FFT
        IQ_fft = np.fft.fft(IQ)
        # Compute power
        power = (np.absolute(IQ_fft) ** 2) / (float(len(IQ)) ** 2)
        # Find frequency of max power
        f_max = f[np.argmax(power)]

        ### refine with continuous FT near first peak
        # frequencies within hwid Hz of peak
        freq = np.arange(f_max-hwid,f_max+hwid+df,df)
        # get continuous FT within hwid Hz of peak
        IQ_fft = np.array([np.sum(IQ*np.exp(-2j*np.pi*time*f)) for f in freq])
        # Compute power
        power = (np.absolute(IQ_fft) ** 2) / (float(len(IQ)) ** 2)
        # Find frequency of max power
        #if abs(f_max - freq[np.argmax(power)]) < hwid:
        f_max = freq[np.argmax(power)]

        return f_max
"""
Revisions
"""
