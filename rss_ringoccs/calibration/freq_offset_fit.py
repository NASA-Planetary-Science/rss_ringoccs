"""
Purpose:
        Compute a fit to the frequency offset residual using
        the frequency offset, predicted sky frequency, and
        reconstructed sky
        frequency.
"""

import numpy as np
import warnings
warnings.simplefilter('ignore', np.RankWarning)
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
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
    Obtains :math:`f(t)_{offset}` from ``calc_freq_offset``,
    :math:`f(t)_{dr}` from ``calc_f_sky_recon``, and :math:`f(t)_{dp}`
    from ``get_f_sky_pred``. Computes a polynomial fit
    :math:`\hat{f}(t)_{resid}` of F-test specified fof_order to
    sigma-clipped residual difference :math:`f(t)_{resid}` between
    observed and predicted frequency offset where the residual is
    given by

    .. math::
        f(t)_{resid} = f(t)_{offset} - (f(t)_{dr}-f(t)_{dp})

    Final frequency offset :math:`\hat{f}(t)_{offset}` is found using
    the polynomial fit :math:`\hat{f}(t)_{resid}` to the frequency
    offset residuals such that

    .. math::
        \hat{f}(t)_{offset} = \hat{f}(t)_{resid} +
        (f(t)_{dr} - f(t)_{dp})

    Arguments:
        :rsr_inst (*object*): object instance of the RSRReader class
        :geo_inst (*object*): object instance of the Geometry class

    Keyword Arguments:
        :f_uso_x (*float*): frequency in Hz of the X-band ultra-stable
                        oscilator onboard the Cassini spacecraft.
                        Default is 8427222034.3405 Hz.
        :verbose (*bool*): when True, enables verbose output mode

    Attributes:
        :f_offset_fit (*np.ndarray*): final frequency offset evaluated using
                        fit to offset residuals
                        :math:`\hat{f}(t)_{offset} =
                        \hat{f}(t)_{resid} + (f(t)_{dr} - f(t)_{dp})`
        :f_spm (*np.ndarray*): SPM at which the offset frequency was sampled
        :f_sky_pred (*np.ndarray*): predicted sky frequency :math:`f(t)_{dp}`
        :f_offset_fit (*np.ndarray*): fit to the residual frequency offset
                        math:`\hat{f}(t)_{resid}` evaluated at
                        ``f_spm``
        :chi_squared (*float*): sum of the squared residual frequency offset fit
                        normalized by the fit value (Pearson's
                        :math:`\chi^2`) such that
                        :math:`\chi^2 = \\frac{1}{N-m}
                        \sum((\hat{f}(t)_{resid}-f(t)_{resid})
                        /\hat{f}(t)_{resid})^2`
                        for :math:`N` data and :math:`m` free
                        parameters (i.e., the polynomial order plus
                        one).
    """

    def __init__(self, rsr_inst, geo_inst, f_uso_x=8427222034.34050,
            verbose=False, write_file=False, fit_resid=False):


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
        #   NOTE: predicted sky frequency extracted from rsr_inst later
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
        inds = [(self.raw_rho>6.25e4)&(self.raw_rho<2.5e5)]
        occ_inds = [(self.raw_rho>7e4)&(self.raw_rho<1.4e5)]

        # Find the max and min SPM values with True boolean indices
        if len(self.raw_spm_vals[inds]) > 2 :
            spm_min = np.min(self.raw_spm_vals[inds])
            spm_max = np.max(self.raw_spm_vals[inds])
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

        # Compute predicted sky frequency
        if verbose:
            print('\tCalculating predicted frequency offset...')
        #spm0, f_sky_pred = rsr_inst.get_f_sky_pred(f_spm=f_spm)

        f_sky_pred = calc_f_sky_recon(f_spm, rsr_inst, sc_name, f_uso,
                kernels)

        # simply fit the offset frequency
        if verbose:
            print('\tCreating sigma clipping mask array...')
        self.__fsr_mask = self.create_mask(f_spm, f_rho, f_offset)

        # if the freq offset covers only a small range (i.e., less than
        # 68% of the data set), then use a low polynomial order to
        # prevent over-fitting of the offset frequency
        if ((f_spm[self.__fsr_mask][-1]-f_spm[self.__fsr_mask][0]) /
            (f_spm[-1]-f_spm[0])) < 0.68 :
            if verbose:
                print('\tInsufficient coverage of occultation. Setting'+
                        '\n\toffset frequency polynomial order to 2.')
            self.poly_order = 2
        # otherwise, compute best polynomial order using F test
        else:
            if verbose:
                print('\tEstimating the best polynomial order...')
            self.poly_order = self.calc_poly_order(f_spm[self.__fsr_mask],
                        f_offset[self.__fsr_mask], verbose=False)

        # Fit frequency offset
        if verbose:
            print('\tCalculating fit to frequency offset...')
        f_offset_fit,chi2 = self.fit_freq_offset(f_spm, f_rho, f_offset)

        # Draw and save reference plot
        if write_file:
            self.plotFORFit(f_spm,f_offset,f_offset_fit,self.__fsr_mask,
                            spm_min,spm_max,geo_inst.t_oet_spm_vals[0],
                            geo_inst.t_oet_spm_vals[-1])

        # set attributes
        self.f_offset_fit = f_offset_fit
        self.f_spm = foff_inst.f_spm
        self.f_offset = foff_inst.f_offset
        self.f_sky_pred  = f_sky_pred
        self.chi_squared = chi2

    def create_mask(self, f_spm, f_rho, f_offset):
        """
        Purpose:
            Creates a Boolean mask array which excludes data based on
            the following critera:
                #. ring or planetary occultation in region prevents
                accurate estimation of the offset frequency
                #. offset frequencies fall more than 5-sigma beyond
                the median offset frequency
                #. adjacent data all excluded by previous requirements
                (excludes noise which by happenstance satisfies the
                above criteria)

        Arguments:
            :f_spm (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        when calculating the offset frequencies for
                        the occultation
            :f_rho (*np.ndarray*): ring intercept radius of the
                        spacecraft signal resampled to match f_spm
            :f_offset (*np.ndarray*): residual sky frequency

        Returns:
            :fsr_mask (*np.ndarray*): Array of booleons, with True for
                                      reliable residual frequency offset.
        """
        dt_spm = round((self.raw_spm_vals[-1]-self.raw_spm_vals[0])
                    /float(len(self.raw_spm_vals)),8)
        dt_fof = round((f_spm[-1]-f_spm[0])/float(len(f_spm)),6)
        df = 12.5*dt_spm*dt_fof

        # Create mask array that includes everything
        fsr_mask = np.array([True for i in range(len(f_offset))],dtype=bool)
        # Compute median, standard deviation, and implememt sigma-clipping
        #   for data which fall in acceptable regions
        fsr_median = np.nanmedian(f_offset[fsr_mask])
        fsr_stdev = 5.*np.sqrt(np.nanmedian(np.square(f_offset-fsr_median)))

        # iteratively check to see if each residual value is within 3 sigma
        for i in range(len(f_offset)):
            # exclude nans
            if np.isnan(f_offset[i]):
                fsr_mask[i] = False
            # exclude values well outside of the possible range
            elif f_offset[i] < -150 or f_offset[i] > 150 :
                fsr_mask[i] = False
            # exclude data outside 3-sigma of median
            elif (f_offset[i] < fsr_median - fsr_stdev) or (
                    f_offset[i] > fsr_median + fsr_stdev):
                fsr_mask[i] = False
            # exclude data with discontinuous jumps
            elif i > 1 and i < len(f_offset)-2:
                chk1 = (abs(f_offset[i]-f_offset[i-1]) > df)
                chk2 = (abs(f_offset[i]-f_offset[i+1]) > df)
                if chk1 or chk2 :
                    fsr_mask[i] = False
                # if this offset freq passes sigma and local variability checks
                # proceed to check variability within the masked data
                if fsr_mask[i] :
                    # look for next reliable offset freq
                    j = int(1)
                    while j < len(f_offset)-(i+j) and not fsr_mask[i+j]:
                        j += 1
                    # look for next reliable offset freq
                    k = int(1)
                    while k < len(f_offset)-(i+j+k) and not fsr_mask[i+j+k]:
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
                            fsr_mask[i+j] = False
            else:
                fsr_mask[i] = False

        # iteratively check adjacent values for false positives
        #   i.e., all four adjacent mask array values are False
        #   first forwards
        for i in range(2,len(fsr_mask)-2):
            if fsr_mask[i]:
                if not fsr_mask[i-2] and not fsr_mask[i-1]:
                    if not fsr_mask[i+1] and not fsr_mask[i+2]:
                        fsr_mask[i] = False

        # now check backwards, just in case false positives were supporting
        # each other and preventing removal
        for i in range(len(fsr_mask)-3,2,-1):
            if fsr_mask[i]:
                if not fsr_mask[i-2] and not fsr_mask[i-1]:
                    if not fsr_mask[i+1] and not fsr_mask[i+2]:
                        fsr_mask[i] = False
                if not fsr_mask[i-1] and not fsr_mask[i+1]:
                    fsr_mask[i] = False

        # if there are no True values in mask array, then reset
        #   and hope for the best
        if not np.any(fsr_mask):
            fsr_mask = np.array([True for i in range(len(f_offset))],
                    dtype=bool)
            for i in range(len(f_offset)):
                # exclude nans
                if np.isnan(f_offset[i]):
                    fsr_mask[i] = False
                # exclude values well outside of the possible range
                elif f_offset[i] < -150 or f_offset[i] > 150 :
                    fsr_mask[i] = False

        # iteratively do polynomial clipping
        for ipol in range(3):
            # compute best polynomial order
            self.poly_order = self.calc_poly_order(f_spm[fsr_mask], f_offset[fsr_mask])
            # Polynomial fit clipping
            pinit = np.polyfit(f_spm[fsr_mask], f_offset[fsr_mask], self.poly_order)

            # Compute standard deviation from fit and implememt sigma-clipping
            #   for data which fall in acceptable regions
            fit_stdev = 5.*np.sqrt(np.nanmedian(np.square(f_offset-
                np.polyval(pinit,f_spm))))

            # if the fit can give us a reasonable constraint, use it to
            #   help sigma clip
            if fit_stdev < 1 :
                # store old mask
                old_mask = fsr_mask
                # Create new mask array that includes everything
                fsr_mask = np.array([True for i in range(len(f_offset))],
                        dtype=bool)
                # iteratively check to see if each residual value is
                # within 5 sigma of fit
                for i in range(len(f_offset)):
                    # exclude values 5-sigma outside of polynomial approximation
                    if (f_offset[i] > np.polyval(pinit,f_spm[i]) -
                            fit_stdev) and (f_offset[i] < np.polyval(
                                pinit,f_spm[i]) + fit_stdev):
                        # exclude nans
                        if np.isnan(f_offset[i]):
                            fsr_mask[i] = False
                        # exclude values well outside of the possible range
                        elif f_offset[i] < -150 or f_offset[i] > 150 :
                            fsr_mask[i] = False
                        # exclude data outside 3-sigma of median
                        elif (f_offset[i] < np.polyval(pinit,f_spm[i]) -
                                fit_stdev) or (f_offset[i] > np.polyval(
                                    pinit,f_spm[i]) + fit_stdev):
                            fsr_mask[i] = False
                        # exclude data with discontinuous jumps
                        elif i > 1 and i < len(f_offset)-2:
                            chk1 = (abs(f_offset[i]-f_offset[i-1]) > df)
                            chk2 = (abs(f_offset[i]-f_offset[i+1]) > df)
                            if chk1 or chk2 :
                                fsr_mask[i] = False
                            # if this offset freq passes sigma and local variability checks
                            # proceed to check variability within the masked data
                            if fsr_mask[i] :
                                # look for next reliable offset freq
                                j = int(1)
                                while j < len(f_offset)-(i+j) and not fsr_mask[i+j]:
                                    j += 1
                                # look for next reliable offset freq
                                k = int(1)
                                while k < len(f_offset)-(i+j+k) and not fsr_mask[i+j+k]:
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
                                        fsr_mask[i+j] = False
                        else:
                            fsr_mask[i] = False
                    else:
                        fsr_mask[i] = False

            # iteratively check adjacent values for false positives
            #   i.e., all four adjacent mask array values are False
            #   first forwards
            for i in range(2,len(fsr_mask)-2):
                if fsr_mask[i]:
                    if not fsr_mask[i-2] and not fsr_mask[i-1]:
                        if not fsr_mask[i+1] and not fsr_mask[i+2]:
                            fsr_mask[i] = False

            # now check backwards, just in case false positives were supporting
            # each other and preventing removal
            for i in range(len(fsr_mask)-3,2,-1):
                if fsr_mask[i]:
                    if not fsr_mask[i-2] and not fsr_mask[i-1]:
                        if not fsr_mask[i+1] and not fsr_mask[i+2]:
                            fsr_mask[i] = False
                    if not fsr_mask[i-1] and not fsr_mask[i+1]:
                        fsr_mask[i] = False

        ## return frequency sky residual mask array
        return fsr_mask

    def calc_poly_order(self, f_spm_cl, f_offset_cl, verbose=False):
        """
        Purpose:
            Use a variant of the F-test to determine the best order
            polynomial to use to fit the frequency offset.

        Arguments:
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
        # dummy F test value
        F = 100.
        if verbose :
            row = [' ORDER I',' SUM SQ RESID',' ORDER I+1',' SUM SQ RESID',' F TEST']
            print('-'*72)
            print("|".join(str(val).ljust(14) for val in row))
            print('-'*72)
        # while F test finds new terms to be significant and the order
        # is at most a ninth order polynomial
        while F > 1. and polyord2 < 11 :
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


    def fit_freq_offset(self, f_spm, f_rho, f_offset,verbose=False):
        """
        Fit a polynomial to residual frequency.

        Arguments:
            :f_spm (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        when calculating the offset frequencies for
                        the occultation
            :f_rho (*np.ndarray*): ring intercept radius of the
                        spacecraft signal resampled to match f_spm
            :f_offset (*np.ndarray*): carrier frequency offset from
                        center of band

        Keyword Arguments:
            :verbose (*bool*): If True, print processing steps

        Returns:
            :f_offset_fit (*np.ndarray*): fit to the residual frequency
                            offset math:`\hat{f}(t)_{resid}` evaluated at
                            ``f_spm``
            :chi2 (*float*): sum of the squared residual frequency offset fit
                            normalized by the fit value (Pearson's
                            :math:`\chi^2`) such that
                            :math:`\chi^2 = \\frac{1}{N-m}
                            \sum((\hat{f}(t)_{resid}-f(t)_{resid})
                            /\hat{f}(t)_{resid})^2`
                            for :math:`N` data and :math:`m` free
                            parameters (i.e., the polynomial order plus
                            one).
        """

        npts = len(f_spm)
        spm_temp = ((f_spm - f_spm[int(npts / 2)])
            / max(f_spm - f_spm[int(npts / 2)]))

        ## fit using polynomial of user-selected order
        coef = np.polyfit(spm_temp[self.__fsr_mask],f_offset[self.__fsr_mask],
                                self.poly_order)

        '''if verbose:
            print('\tPolynomial sum squared residuals:',stats[0])'''

        f_offset_fit = np.polyval( coef, spm_temp )
        v = float(len(f_offset[self.__fsr_mask])) - (self.poly_order+1)
        chi2 = np.sum(np.square(f_offset_fit[self.__fsr_mask]-
            f_offset[self.__fsr_mask]) / f_offset_fit[self.__fsr_mask])

        return f_offset_fit,chi2

    # Create and save a plot of the offset residual fit
    def plotFORFit(self,spm,resid,fit,mask,spm_min,spm_max,occ_min,occ_max):
        """
        Purpose:
            Plot results of the automated frequency offset residual
            fit and save plot to a file. File name will match the
            .LBL and .TAB nomenclature.

        Arguments:
            :spm (*np.ndarray*): SPM sampled by ``calc_freq_offset``
                        when calculating the offset frequencies for
                        the occultation
            :resid (*np.ndarray*): residual sky frequency
            :fit (*np.ndarray*): polynomial fit to the residual sky
                        frequency
            :mask (*np.ndarray*): boolean array used to mask residual
                        sky frequency for the polynomial fitting
            :spm_min (*float*): start of occultation in SPM
            :spm_max (*float*): end of occultation in SPM
        """
        #generate plot file names
        filenames,outdirs = construct_filepath(self.rev_info,'FORFIT')
        # set up subplot
        ax = plt.figure().add_subplot(111)
        # residuals used for fit
        plt.plot(spm[mask],resid[mask],'.k')
        # all residuals
        plt.plot(spm,resid,'-',color='0.5',lw=1)
        # indicate limits for ring system
        #plt.axvline(occ_min,dashes=[12,4],color='0.2')
        #plt.axvline(occ_max,dashes=[12,4],color='0.2')
        # fit to residuals
        plt.plot(spm,fit,'-r')
        # limits to plot
        plt.xlim(spm_min-100,spm_max+100)
        plt.ylim(np.nanmin(resid[mask])-0.1,np.nanmax(resid[mask])+0.1)
        # labels
        plt.xlabel('SPM (sec)')
        plt.ylabel(r'$f_{predict}-f_{observe}$')
        plt.text(0.4,0.95,'PolyOrder: '+str(self.poly_order),transform =
                ax.transAxes)
        # output
        for file,dir in zip(filenames,outdirs):
            plt.title(file)
            outfile = dir + file + '.PDF'
            print('\tSaving frequency offset fit plot to: \n\t\t' + '/'.join(outfile.split('/')[0:5]) + '/\n\t\t\t' + '/'.join(outfile.split('/')[5:]))
            plt.savefig(outfile)
        plt.close()
"""
Revisions:
"""
