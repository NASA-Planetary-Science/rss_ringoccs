"""
Purpose:
        Compute a fit to the frequency offset using offset frequencies
        calculated from raw data, sigma-clipping frequencies
        contaminated by rings, and fitting with a polynomial of order
        determined by an iterative F-test.
"""

import numpy as np
import warnings
warnings.simplefilter('ignore', np.RankWarning)
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.signal import savgol_filter,argrelmax
import sys

from .calc_f_sky_recon import calc_f_sky_recon
from .calc_freq_offset import calc_freq_offset
from ..tools.write_output_files import construct_filepath

import sys
sys.path.append('../../')
import rss_ringoccs as rss
sys.path.remove('../../')


class FreqOffsetFit(object):
    """
    :Purpose:
        Obtains :math:`f(t)_{offset}` from ``calc_freq_offset``,
        :math:`f(t)_{dr}` from ``calc_f_sky_recon``. Computes a polynomial
        fit :math:`\hat{f}(t)_{offset}` of F-test specified order to
        sigma-clipped frequency offset. Final sky frequency
        :math:`\hat{f}(t)_{sky}` is calculated by summing the polynomial
        fit :math:`\hat{f}(t)_{offset}` with the reconstructed sky
        frequency :math:`f(t)_{dr}`.
    :Arguments:
        :rsr_inst (*object*): object instance of the RSRReader class
        :geo_inst (*object*): object instance of the Geometry class
    :Keyword Arguments:
        :f_uso_x (*float*): frequency in Hz of the X-band ultra-stable
                        oscilator onboard the Cassini spacecraft.
                        Default is 8427222034.3405 Hz.
        :verbose (*bool*): when True, enables verbose output mode
    :Attributes:
        :f_offset_fit (*np.ndarray*): fit to frequency offset :math:`\hat{f}(t)_{offset}
        :f_spm (*np.ndarray*): SPM at which the offset frequency was sampled
        :f_sky_recon (*np.ndarray*): reconstructed sky frequency :math:`f(t)_{dr}`
        :f_offset_fit (*np.ndarray*): fit to the frequency offset
                        math:`\hat{f}(t)_{offset}` evaluated at ``f_spm``
        :chi_squared (*float*): sum of the squared residual difference between
                        the frequency offset and the frequency offset fit
                        normalized by the fit value (Pearson's
                        :math:`\chi^2`) such that
                        :math:`\chi^2 = \\frac{1}{N-m}
                        \sum((\hat{f}(t)_{offset}-f(t)_{offset})
                        /\hat{f}(t)_{offset})^2`
                        for :math:`N` data and :math:`m` free
                        parameters (i.e., the polynomial order plus
                        one).
    """

    def __init__(self, rsr_inst, geo_inst, f_uso_x=8427222034.34050,
            verbose=False, write_file=False, fof_lims=None):


        # Check inputs for validity
        if not isinstance(rsr_inst, rss.rsr_reader.RSRReader):
            sys.exit('ERROR (FreqOffsetFit): rsr_inst input must be an '
                + 'instance of the RSRReader class')

        if not isinstance(geo_inst, rss.occgeo.Geometry):
            sys.exit('ERROR (FreqOffsetFIt): geo_inst input must be an '
                + 'instance of the Geometry class')

        if not isinstance(verbose, bool):
            print('WARNING (FreqOffsetFit): verbose input must be boolean. '
                + 'Ignoring current input and setting to False')
            verbose = False

        # Extract necessary information from input instances
        self.band = rsr_inst.band
        self.year = rsr_inst.year
        self.doy = rsr_inst.doy
        self.dsn = rsr_inst.dsn
        self.raw_spm_vals = rsr_inst.spm_vals
        self.__IQ_m = rsr_inst.IQ_m
        self.rev_info = rsr_inst.rev_info


        spm_geo = geo_inst.t_oet_spm_vals
        rho_geo = geo_inst.rho_km_vals
        kernels = geo_inst.kernels
        self.profdir = geo_inst.get_profile_dir()
        self.rev_info = geo_inst.rev_info
        sc_name = geo_inst.history['Positional Args']['spacecraft']

        # Adjust USO frequency by wavelength
        if self.band == 'X':
            f_uso = f_uso_x
        elif self.band == 'S':
            f_uso = f_uso_x*(3.0/11.0)
        elif self.band == 'K':
            f_uso = f_uso_x*3.8
        else:
            raise ValueError('WARNING (freq_offset_fit.py): Invalid frequency '
                            + 'band!')

        # Compute spline coefficients relating SPM to rho
        rho_geo_spl_coef = splrep(spm_geo, rho_geo)

        # compute max and min SPM values for occultation
        # Evaluate spm-to-rho spline at raw SPM to get raw rho sampling
        #    that matches SPM values
        self.raw_rho = splev(self.raw_spm_vals,rho_geo_spl_coef)

        # Create boolean mask where True is within occultation range and
        #    False is outside the occultation range -- this is generalized
        #    to work for diametric and chord occultations
        inds = [(self.raw_rho>6.5e4)&(self.raw_rho<1.75e5)]
        #occ_inds = [(self.raw_rho>7e4)&(self.raw_rho<1.4e5)]

        # If spm limits are provided in the right format, use those
        # to select a portion of the occultation
        if fof_lims != None and len(fof_lims) == 2 and fof_lims[1]-fof_lims[0] > 100.:
            spm_min = fof_lims[0]
            spm_max = fof_lims[1]
        # Find the max and min SPM values with True boolean indices
        else:
            # account for Rev 58 being awful -- chord lacking all but C ring in egress
            if self.year == 2008 and self.doy == 39:
                spm_min = 6.4e4
                spm_max = 6.7e4
            # if radius range includes sufficient data
            elif len(self.raw_spm_vals[inds]) > 2 :
                spm_min = np.min(self.raw_spm_vals[inds])
                spm_max = np.max(self.raw_spm_vals[inds])
            # if insufficient data in radial range, use full SPM range
            else:
                print('Error in estimating SPM range for frequency offset!')
                spm_min = self.raw_spm_vals[0]
                spm_max = self.raw_spm_vals[-1]

        # Calculate offset frequency within given SPM limits
        if verbose:
            print('\tCalculating observed frequency offset...')

        foff_inst = calc_freq_offset(rsr_inst,spm_min,spm_max)
        f_spm, f_offset = foff_inst.f_spm, foff_inst.f_offset

        # Interpolate rho to frequency time values
        f_rho = splev(f_spm, rho_geo_spl_coef)

        # Compute reconstructed sky frequency
        if verbose:
            print('\tCalculating reconstructed frequency offset...')

        f_sky_recon = calc_f_sky_recon(f_spm, rsr_inst, sc_name, f_uso,
                kernels)

        # simply fit the offset frequency
        if verbose:
            print('\tCreating sigma clipping mask array...')
        self.__mask = self.create_mask(f_spm, f_rho, f_offset)

        # Fit frequency offset
        if verbose:
            print('\tCalculating fit to frequency offset...')
        f_offset_fit,chi2 = self.fit_freq_offset(f_spm, f_rho, f_offset)

        # Draw and save reference plot
        if write_file:
            self.plotFORFit(f_spm,f_offset,f_offset_fit,self.__mask,
                            spm_min,spm_max,geo_inst.t_oet_spm_vals[0],
                            geo_inst.t_oet_spm_vals[-1])

        # set attributes
        self.f_offset_fit = f_offset_fit
        self.f_spm = foff_inst.f_spm
        self.f_offset = foff_inst.f_offset
        self.f_sky_recon  = f_sky_recon
        self.chi_squared = chi2

    def create_mask(self, f_spm, f_rho, f_offset, polyclip=False, Cval=7.5):#12.5):
        """
        :Purpose:
            Creates a Boolean mask array which excludes data based on
            the following critera:
                #. ring or planetary occultation in region prevents
                   accurate estimation of the offset frequency
                #. offset frequencies fall more than 5-sigma beyond
                   the median offset frequency
                #. offset frequencies vary by more than  0.25 Hz relative
                   to neighboring  offset frequencies
                #. adjacent data all excluded by previous requirements
                   (excludes noise which by happenstance satisfies the
                   above criteria)
        :Arguments:
            :f_spm (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        when calculating the offset frequencies for
                        the occultation
            :f_rho (*np.ndarray*): ring intercept radius of the
                        spacecraft signal resampled to match f_spm
            :f_offset (*np.ndarray*): frequency offset
        :Keyword Arguments:
            :Cval (*float*): constant scale factor which sets the
                        tolerance threshold for outliers relative
                        to neighboring frequency offsets. Scale factor
                        for 1 and 16 kHz files is such that
                        the frequency threshold is 0.02*Cval in Hz.
                        Default is 5.
            :polyclip (*bool*): boolean specifying whether to do a
                        3rd-order polynomial fit to the sigma-clipped
                        frequency offsets and perform additional
                        sigma-clipping based on the results of the fit.
                        Default is False. It is highly recommended
                        that users keep the default unless the final
                        frequency offset fit results in extremely poor
                        phase drift corrections.
        :Returns:
            :mask (*np.ndarray*): Array of booleons, with True for
                                      reliable frequency offset.
        """
        dt_spm = round((self.raw_spm_vals[-1]-self.raw_spm_vals[0])
                    /float(len(self.raw_spm_vals)),8)
        dt_fof = round((f_spm[-1]-f_spm[0])/float(len(f_spm)),6)
        if dt_spm < 1e-3 :
            df = Cval*dt_fof*1e-3
        else:
            df = Cval*dt_spm*dt_fof

        # Compute median, standard deviation, and implememt sigma-clipping
        #   for data which fall in acceptable regions
        f_median = np.nanmedian(f_offset)
        f_stdev = 5.*np.sqrt(np.nanmedian(np.square(f_offset-f_median)))
        if f_stdev > 15.:
            f_stdev = 15.
        # difference between offset frequency and its median
        median_diff = abs(f_offset-f_median)

        # iteratively check to see if each freq offset value is within
        # specified sigma (here, 5 standard deviations)
        mask = self.__sigma_clip(f_spm,f_offset,median_diff,f_stdev,df=df)
        # check for false positives
        mask = self.__neighbor_check(mask)

        # if there are no True values in mask array, then exclude nans
        # and hope for the best
        if len(f_offset[mask]) < 5:
            mask = np.array([True]*len(f_offset),dtype=bool)
            for i in range(len(f_offset)):
                # exclude nans
                if np.isnan(f_offset[i]):
                    mask[i] = False
                # exclude values well outside of the possible range
                elif f_offset[i] < -150 or f_offset[i] > 150 :
                    mask[i] = False
                # exclude B ring
                elif f_rho[i] > 9.2e4 and f_rho[i] < 1.18e5 :
                    mask[i] = False
            # check for false positives
            mask = self.__neighbor_check(mask)

        # iteratively do polynomial clipping
        if polyclip:
            for ipol in range(3):
                # polynomial fit
                pinit = np.polyfit(f_spm[mask], f_offset[mask], 2)

                # Compute standard deviation from fit and implememt sigma-clipping
                #   for data which fall in acceptable regions
                fit_stdev = 5.*np.sqrt(np.nanmedian(np.square(f_offset-
                    np.polyval(pinit,f_spm))))

                # if the fit can give us a reasonable constraint, use it to
                #   help sigma clip
                if fit_stdev < 1 :
                    # store old mask
                    old_mask = mask
                    # get absolute difference between polynomial and data
                    poly_diff = abs( f_offset - np.polyval(pinit,f_spm) )
                    # iteratively check to see if each freq offset value is
                    # within 5 sigma of fit
                    new_mask = self.__sigma_clip(f_offset,poly_diff,fit_stdev,df=df)

                    # check for false positives
                    mask = self.__neighbor_check(new_mask)

        ## return frequency offset mask array
        return mask

    def __sigma_clip(self, f_spm, f_offset, diff, stdev, df=0.25 ):
        # starting mask that is all-inclusive
        mask = np.array([True]*len(f_offset),dtype=bool)
        # iteratively check to see if each freq offset value is
        # within appropriate distance of fit
        for i in range(len(f_offset)):
            # exclude nans
            if np.isnan(f_offset[i]):
                mask[i] = False
            # exclude values well outside of the possible frequency range
            elif f_offset[i] < -150 or f_offset[i] > 150 :
                mask[i] = False
            # exclude data more than 5 sigma outside projected trend
            elif diff[i] > stdev :
                mask[i] = False
            # exclude data with variability greater than df
            elif i > 1 and i < len(f_offset)-2 :
                chk1 = (abs(f_offset[i]-f_offset[i-1]) > df)
                chk2 = (abs(f_offset[i]-f_offset[i+1]) > df)
                if chk1 or chk2 :
                    mask[i] = False
                # if this offset freq passes sigma and local variability checks
                # proceed to check variability within the masked data
                else :
                    # look for next reliable offset freq
                    j = int(1)
                    while j < len(f_offset)-(i+j) and not mask[i+j]:
                        j += 1
                    # look for next reliable offset freq
                    k = int(1)
                    while k < len(f_offset)-(i+j+k) and not mask[i+j+k]:
                        k += 1
                    # compare difference
                    chk3 = (abs(f_offset[i]-f_offset[i+j]) > df)
                    # compare difference
                    chk4 = (abs(f_offset[i+j]-f_offset[i+j+k]) > df)
                    # and have been defined
                    if (i+j) < (i+j+k) and (i+j) < len(f_offset) and (i+j+k) < len(f_offset) :
                        # if both checks pass
                        if chk3 and chk4 :
                            # data is an outlier, so exclude
                            mask[i+j] = False

        return mask

    def __neighbor_check(self, mask):
        # iteratively check adjacent values for false positives
        #   i.e., all four adjacent mask array values are False
        #   first forwards
        for i in range(2,len(mask)-2):
            if mask[i]:
                if not mask[i-2] and not mask[i-1]:
                    if not mask[i+1] and not mask[i+2]:
                        mask[i] = False
        # now check backwards, just in case false positives were supporting
        # each other and preventing removal
        for i in range(len(mask)-3,2,-1):
            if mask[i]:
                if not mask[i-2] and not mask[i-1]:
                    if not mask[i+1] and not mask[i+2]:
                        mask[i] = False
                if not mask[i-1] and not mask[i+1]:
                    mask[i] = False
        return mask

    def calc_poly_order(self, f_spm_cl, f_offset_cl, verbose=False, max_order=9):
        """
        Use a variant of the F-test to determine the best order
        polynomial to use to fit the frequency offset.
        Arguments
            :f_spm_cl (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        and clipped by the initial boolean mask.
            :f_offset_cl (*np.ndarray*): carrier frequency offset from
                        center of band
        """
        # seed polynomial orders
        polyord1 = 0
        polyord2 = 1
        # length of data set
        n = float(len(f_offset_cl))
        # make sure polynomial order doesn't exceed length of data set
        if n < max_order :
            max_order = n-1
        if max_order < 1 :
            max_order = 1
        # dummy F test value
        F = 100.
        if verbose :
            row = [' ORDER I',' SUM SQ RESID',' ORDER I+1',' SUM SQ RESID',' F TEST']
            print('-'*72)
            print("|".join(str(val).ljust(14) for val in row))
            print('-'*72)
        # while F test finds new terms to be significant and the order
        # is at most a ninth order polynomial
        while F > 1. and polyord2 < max_order+1 :
            # if new term is significant update and continue
            polyord1 = polyord2
            polyord2 += 1
            # compute polynomial orders
            p1 = np.polyfit(f_spm_cl,f_offset_cl,polyord1)
            p2 = np.polyfit(f_spm_cl,f_offset_cl,polyord2)
            # get degrees of freedom
            df1 = n-float(len(p1))
            df2 = n-float(len(p2))
            # compute fitted values
            y1 = np.polyval(p1,f_spm_cl)
            y2 = np.polyval(p2,f_spm_cl)
            # compute sum of squared differences
            sqsum1 = np.sum(np.square(y1 - f_offset_cl))
            sqsum2 = np.sum(np.square(y2 - f_offset_cl))
            # F test
            num = (sqsum1-sqsum2)/(df1-df2)
            den = sqsum2/df2
            F = num/den
            if verbose :
                row = [polyord1,round(sqsum1,4),polyord2,round(sqsum2,4),round(F,4)]
                print("|".join(' '+str(val).ljust(13) for val in row))
        if verbose :
            print('-'*72)
        # pass back the final polynomial order
        return polyord1

    def get_disconts(self,f_spm,f_rho,f_offset,df=0.1):
        """
        :Purpose:
            Find any possible discontinuities in the frequency offset by finding
            local maxima in the second derivative of the masked (sigma-clipped)
            frequency offset data. These maxima correspond to sharp edges in the
            frequency offset data. The second derivative is computed using a
            central finite difference approximation. Discontinuities are found
            using a relative maximum algorithm implemented by scipy.signal in
            the method argrelmax to find the peaks in the absolute value of the
            second derivative (these will correspond to suddent changes in slope).
            Discontinuities found by this method are then checked to remove
            false positives, rejecting cases of obstruction by the B ring and
            segments containing less than 10 data.
        :Arguments:
            :f_spm (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        when calculating the offset frequencies for
                        the occultation
            :f_rho (*np.ndarray*): ring intercept radius of the
                        spacecraft signal resampled to match f_spm
            :f_offset (*np.ndarray*): carrier frequency offset from
                        center of band
        """
        # second derivative with central finite difference
        df_cfd = self.__centdiff2(f_spm[self.__mask],f_offset[self.__mask])
        # take absolute value
        #abs_d2fdt2 = abs(df_savgol)
        abs_d2fdt2 = abs(df_cfd)
        # find peaks in 2nd derivative
        # (absolute 2nd derivative will spike at an "instantaneous change", i.e.,
        # where a jump, cusp, or other discontinuity occurs)
        rmax, = argrelmax(abs_d2fdt2)
        # set cutoff in derivative to 2/3 max
        cut = np.max(abs_d2fdt2)/1.5
        # & abs_d2fdt2[d]>1e3
        disconts = [0]+[d for d in rmax if abs_d2fdt2[d]>cut]+[-1]
        '''t_dc = [f_spm[self.__mask][d] for d in disconts]
        df_dc = [abs_d2fdt2[d] for d in disconts]
        #plt.plot(t,abs(df_savgol))
        plt.plot(f_spm[self.__mask],abs(df_cfd))#*np.max(abs(df_savgol))/np.max(abs(df_cfd)))
        plt.plot(t_dc,df_dc,'.r')
        plt.show()'''

        # check for spurious discontinuities
        rmv = [] # array of discontinuities to remove
        for i in range(1,len(disconts)-1):
            # make sure segment is populated by sufficient data
            check1 = len(f_offset[self.__mask][disconts[i-1]:disconts[i]])<10
            check2 = len(f_offset[self.__mask][disconts[i]:disconts[i+1]])<10
            # make sure discontinuity is not a noisy spike
            abdif = [abs(f_offset[self.__mask][disconts[i]]-f_offset[self.__mask][disconts[i]-1]),
                     abs(f_offset[self.__mask][disconts[i]+1]-f_offset[self.__mask][disconts[i]]),
                     abs(f_offset[self.__mask][disconts[i]+1]-f_offset[self.__mask][disconts[i]-1])]
            check3 = (abdif[0]<df)&(abdif[1]<df)&(abdif[2]<df)
            # make sure discontinuity is not a false positive due to the B ring
            check4 = ((f_rho[self.__mask][disconts[i]]>9.1e4)&(f_rho[self.__mask][disconts[i]]<1.19e5))
            #check5 = ((abs(f_rho[self.__mask][disconts[i]+1]-f_rho[self.__mask][disconts[i]])>2e4)|
            #         (abs(f_rho[self.__mask][disconts[i]]-f_rho[self.__mask][disconts[i]-1])>2e4))
            #print(check1,check2,check3,check4,disconts[i],f_spm[self.__mask][disconts[i]],f_rho[self.__mask][disconts[i]])
            # if B ring or insufficient data or false positive
            if check1 or check2 or check3 or check4 :
                # flag for removal
                rmv += [disconts[i]]
        # remove flagged discontinuities from discontinuity list
        for r in rmv:
            disconts.remove(r)
        # store as attribute, including start and finish of data set
        self.disconts = disconts
        return

    # calculate central finite difference as approximation of 2nd order derivative
    def __centdiff2(self,x,f):
        xb = x[:-2]  # x-h, "backwards"
        xi = x[1:-1] # x
        xf = x[2:]   # x+h, "forward"
        fb = f[:-2]  # f(x-h), "backwards"
        fi = f[1:-1] # f(x)
        ff = f[2:]   # f(x+h), "forwards"
        h = (xf-xb)/2      # h as average 0.5*[ ((x+h)-x) + (x-(x-h)) ]
        h2 = (xi-xb)*(xf-xi) # h squared as backward difference * forward difference
        df2dx2 = (ff-2*fi+fb)/h2
        return np.concatenate([[0],df2dx2,[0]])


    def fit_freq_offset(self, f_spm, f_rho, f_offset,verbose=False):
        """
        :Purpose:
            Fit a polynomial to each frequency offset segment determined
            by sigma clipping and searching for discontinuities. Order
            of the polynomial is determined for each segment by an iterative
            F test.
        :Arguments:
            :f_spm (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        when calculating the offset frequencies for
                        the occultation
            :f_rho (*np.ndarray*): ring intercept radius of the
                        spacecraft signal resampled to match f_spm
            :f_offset (*np.ndarray*): carrier frequency offset from
                        center of band
        :Keyword Arguments:
            :verbose (*bool*): If True, print processing steps
        :Returns:
            :f_offset_fit (*np.ndarray*): fit to the frequency
                            offset math:`\hat{f}(t)_{offset}` evaluated at
                            ``f_spm``
            :chi2 (*float*): sum of the squared residual difference between
                            frequency offset and frequency offset fit
                            normalized by the fit value (Pearson's
                            :math:`\chi^2`) such that
                            :math:`\chi^2 = \\frac{1}{N-m}
                            \sum((\hat{f}(t)_{offset}-f(t)_{offset})
                            /\hat{f}(t)_{offset})^2`
                            for :math:`N` data and :math:`m` free
                            parameters (i.e., the polynomial order plus
                            one).
        """
        npts = len(f_spm)
        spm_temp = ((f_spm - f_spm[int(npts / 2)])
            / max(f_spm - f_spm[int(npts / 2)]))
        f_offset_fit = np.zeros(len(f_offset))

        # prior to USO failure, no discontinuities exists
        if self.year < 2011 :
            # infer polynomial order
            poly_order = self.calc_poly_order(spm_temp[self.__mask],f_offset[self.__mask])
            self.poly_order = poly_order
            # fit with polynomial using inferred polynomial
            coef = np.polyfit(spm_temp[self.__mask],f_offset[self.__mask],poly_order)
            # store fit
            f_offset_fit = np.polyval( coef, spm_temp )
            # store some attributes
            self.disconts=[0,-1]
            # fit assessment statistics
            v = float(len(f_offset[self.__mask])) - (np.max(self.poly_order)+1)
            chi2 = np.sum(np.square(f_offset_fit[self.__mask]-
                f_offset[self.__mask]) / f_offset_fit[self.__mask])

        # find discontinuities in post-USO frequency offset before fitting
        else:
            # get discontinuities from second derivative
            self.get_disconts(spm_temp,f_rho,f_offset)
            # store polynomial orders for each segment
            self.poly_order = []
            # fit each segment by iterating over discontinuities
            for i in range(len(self.disconts)-1):
                poly_order = self.calc_poly_order(
                                  spm_temp[self.__mask][self.disconts[i]+1:self.disconts[i+1]],
                                  f_offset[self.__mask][self.disconts[i]+1:self.disconts[i+1]],
                                  max_order=4)
                ## fit using polynomial of user-selected order
                coef = np.polyfit(spm_temp[self.__mask][self.disconts[i]+1:self.disconts[i+1]],
                                  f_offset[self.__mask][self.disconts[i]+1:self.disconts[i+1]],
                                  poly_order)
                ## find spm limits
                if i == 0:
                    # if no discontinuities exist, fit entire freq offset profile
                    if len(self.disconts)==2:
                        spm_min = spm_temp[0]
                        spm_max = spm_temp[-1]
                    # otherwise, store fit starting at initial SPM
                    else:
                        spm_min = spm_temp[0]
                        spm_max = spm_temp[self.__mask][self.disconts[i+1]+1]
                # if at the final segment, store fit ending at final SPM
                elif i == len(self.disconts)-2:
                    spm_min = spm_temp[self.__mask][self.disconts[i]-1]
                    spm_max = spm_temp[-1]
                # otherwise, define segment extent from one discontinuity to the next
                else:
                    spm_min = spm_temp[self.__mask][self.disconts[i]-1]
                    spm_max = spm_temp[self.__mask][self.disconts[i+1]+1]
                ## clip spm to given limits
                spm_clp = [(spm_temp>=spm_min)&(spm_temp<=spm_max)]
                ## evaluate fit over given SPM range
                f_offset_fit[spm_clp] = np.polyval( coef, spm_temp[spm_clp] )
                # store polynomial order for this segment
                self.poly_order += [poly_order]
            # net fit statistics over all segments
            v = float(len(f_offset[self.__mask])) - (np.max(self.poly_order)+1)
            chi2 = np.sum(np.square(f_offset_fit[self.__mask]-
                f_offset[self.__mask]) / f_offset_fit[self.__mask])
        # return fit and fit assessment
        return f_offset_fit,chi2

    # Create and save a plot of the freq offset fit
    def plotFORFit(self,spm,f_offset,fit,mask,spm_min,spm_max,occ_min,occ_max):
        """
        :Purpose:
            Plot results of the automated frequency offset
            fit and save plot to a file. File name will match the
            .LBL and .TAB nomenclature.
        :Arguments:
            :spm (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        when calculating the offset frequencies for
                        the occultation
            :f_offset (*np.ndarray*): frequency offset
            :fit (*np.ndarray*): polynomial fit to the frequency offset
            :mask (*np.ndarray*): boolean array used to mask frequency
                        offset for the polynomial fitting
            :spm_min (*float*): start of occultation in SPM
            :spm_max (*float*): end of occultation in SPM
        """
        #generate plot file names
        filenames,outdirs = construct_filepath(self.rev_info,'FORFIT')
        # set up subplot
        ax = plt.figure().add_subplot(111)
        # frequency offsets used for fit
        plt.plot(spm[mask],f_offset[mask],'.k')
        # all frequency offsets
        plt.plot(spm,f_offset,'-',color='0.5',lw=1)
        # fit to frequency offset
        plt.plot(spm,fit,'-r')
        # discontinuities
        for d in self.disconts:
            plt.axvline(spm[mask][d],color='C2')
        # limits to plot
        plt.xlim(spm_min-100,spm_max+100)
        plt.ylim(np.nanmin(f_offset[mask])-0.1,np.nanmax(f_offset[mask])+0.1)
        # labels
        plt.xlabel('SPM (sec)')
        plt.ylabel(r'Frequency Offset (Hz)')
        bbox = dict(facecolor='white',edgecolor='none',alpha=0.75,boxstyle='round')
        plt.text(0.4,0.95,'PolyOrder: '+str(self.poly_order),transform =
                ax.transAxes,bbox=bbox)
        # output
        for file,dir in zip(filenames,outdirs):
            plt.title(file)
            outfile = dir + file + '.PDF'
            plt.savefig(outfile)
            print('\tFrequency offset fit plot saved to: ' + outfile)
        plt.close()
