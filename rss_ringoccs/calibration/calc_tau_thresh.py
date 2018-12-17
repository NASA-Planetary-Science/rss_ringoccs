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
from scipy.signal import spectrogram,resample_poly
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

        # convert measured IQ to power at the spm_cal sampling
        p_m = np.interp(cal_inst.t_oet_spm_vals,rsr_inst.spm_vals,abs(rsr_inst.IQ_m**2))
        dt = cal_inst.t_oet_spm_vals[1]-cal_inst.t_oet_spm_vals[0]

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
        noise = self.find_noise(rho_km,p_m,cal_inst.p_free_vals,dt)
        # signal
        #signal = np.interp(spm_spec, cal_inst.t_oet_spm_vals, cal_inst.p_free_vals)
        signal = cal_inst.p_free_vals
        #
        bandwidth = abs( rho_dot_kms / res_km )
        snr = signal/noise

        tau = -np.sin(B_rad)*np.log(0.5*Calpha*bandwidth*snr)

        # compute SNR and set attribute
        self.snr = snr

        ###
        ### ~~ COMPUTE THRESHOLD OPTICAL DEPTH ~~
        ###
        self.tau_thresh = tau

    def find_noise(self,rho,p_m,p_fit,dt):
        """
        Purpose:
            Locate the additive receiver noise within the data set
            based on initial estimate that SNR > 10 and the
            assumption that the receiver noise will not appear within
            the radius range of the ring system in the occultation.

        Arguments:
            :rho (*np.ndarray*): Radius sampled at the calibration
                        SPM sampling rate (``dt_cal`` in the
                        calibration object class)
            :p_m (*np.ndarray*): complex measured signal power
            :p_fit (*np.ndarray*): polynomial fit to the freespace
                        power; representative of the expected
                        intrinsic spacecraft power
            :dt (*float*): sampling rate of the SPM (``dt_cal``) for
                        use in determining the indices of the first
                        and last 1,000 seconds of the occultation in
                        order to determine the location of receiver
                        noise within the data set
        """
        # snr threshold
        thresh = 0.1 * np.nanmin(p_fit[(rho>7e4)&(rho<1.4e5)])
        # storage list
        p_noise = []
        # for power values in first 1,000 secs of data set
        for i in range(1,int(1e3/dt)):
            # if in the right SNR ballpark
            if p_m[i-1] < thresh and p_m[i] < thresh and p_m[i+1] < thresh :
                # store for statistics
                p_noise += [p_m[i]]
        # for power values in last 1,000 secs of data set
        for i in range(len(p_m)-int(1e3/dt),len(p_m)-1):
            # if in the right SNR ballpark
            if p_m[i-1] < thresh and p_m[i] < thresh and p_m[i+1] < thresh :
                # store for statistics
                p_noise += [p_m[i]]
        # return median of the noise power
        return np.nanmedian(p_noise)

"""
Revision history
    Nov 29 2018 - sflury -- original
"""
