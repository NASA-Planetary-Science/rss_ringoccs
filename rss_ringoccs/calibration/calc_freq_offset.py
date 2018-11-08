import numpy as np

'''
    calc_freq_offset

        Class for computing the frequency corresponding to the maximum power in
        the continuous FFT power spectrum

    Nov 01 2018 - sflury        -- original, used components of v1.0 script
'''
class calc_freq_offset(object):
    '''
        __init__

            Takes
                self
                rsr_inst

            Optional
                dt_freq         -- half the width of the FFT window
                delta_t_cent    -- incremental change in SPM window center
                                   Note: this will determine the speed of the
                                   frequency offset calculations. Window spacing
                                   less than 200 seconds will take time:
                                        100 sec spacing      4 sec
                                         50 sec spacing      8 sec
                                         25 sec spacing     16 sec
                                         10 sec spacing     46 sec
                                    Default is 20 sec spacing.


            Calls functions to set up and run continuous FFTs and get peak freqs
    '''
    def __init__(self,rsr_inst,dt_freq=131.072,n_FFTs=5e2):

        # Get raw SPM and raw I & Q from RSR instance
        self.spm_vals = rsr_inst.spm_vals
        self.IQ_m = rsr_inst.IQ_m
        # raw SPM sampling rate in seconds
        self.dt = np.nanmedian([self.spm_vals[i+1] - self.spm_vals[i] for i in range(len(self.spm_vals)-1)])

        # set frequency sampling as attribute
        self.dt_freq = dt_freq
        # set number of samples
        self.n_FFTs = n_FFTs
        # finds and sets ranges/limits for FFTs
        ### ~~ WARNING ~~ ###
        # This will clip the SPM and IQ data attributes
        #self.__set_fft_ranges()

        # finds offset frequencies
        self.__find_offset_freqs()



    '''
        __find_freqs

            Takes
                self

            Iteratively calls __find_peak_freq for slices of SPM and IQself.

            Sets
                f_spm           -- frequency offset files sampled wrt SPM
    '''
    def __find_offset_freqs(self):

        #
        delta_t_cent = (np.max(self.spm_vals)-np.min(self.spm_vals))/float(self.n_FFTs)

        # storage lists -- later converted to arrays and stored as attributes
        spms = []
        freqs = []
        # iteratively compute peak frequency
        for spm_mid in np.arange(np.min(self.spm_vals),np.max(self.spm_vals)+delta_t_cent,delta_t_cent):
            # set indices using boolean mask over a range of SPM
            ind = [(self.spm_vals>=spm_mid)&(self.spm_vals<spm_mid+self.dt_freq)]
            # make sure data are included in range
            if len(self.spm_vals[ind]) > 2 :
                # get average SPM in window
                spms += [np.nanmean(self.spm_vals[ind])]
                # get frequency of peak in power spectrum
                freqs += [self.__find_peak_freq(self.IQ_m[ind])]
        # convert to arrays and store as attributes
        self.f_spm = np.array(spms)
        self.f_offset = np.array(freqs)+0.01
        # pay no attention to the extra 0.01 correction needed but inexplicable

    '''
        __find_peak_freq

            Takes
                spm         -- spm_vals within the current window
                IQ          -- IQ_m vals within the current window

            Computes continuous FFT, finds frequency at max power

            Returns
                f_max       -- frequency at max power
    '''
    def __find_peak_freq(self,IQ):

        # Compute and apply Hamming window
        weight = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(float(len(IQ))) /
            (float(len(IQ)) - 1)))
        IQ *= weight

        # Get FFT frequencies
        f = np.fft.fftfreq(IQ.size)/self.dt
        # Compute FFT with extra exp( -j * 2pi * f * dt ) factor for continuous
        IQ_fft = np.fft.fft(np.exp(-1j*2*np.pi*f*self.dt)*IQ)
        # Compute power
        power = (np.absolute(IQ_fft) ** 2) / (float(len(IQ)) ** 2)

        # Find frequency of max power
        f_max = f[np.argwhere(power==np.nanmax(power))][0][0]

        return f_max
