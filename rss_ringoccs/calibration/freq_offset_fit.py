"""

freq_offset_fit_jol.py

Purpose: Makes a fit to the frequency offset made from freq_offset.py, using
         the frequency offset, predicted sky frequency, reconstructed sky
         frequency, and a fit to residual frequency.
"""

import numpy as np
#from numpy.polynomial import polynomial as poly
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
import sys
import pdb
import pickle

from .calc_f_sky_recon import calc_f_sky_recon
#from .calc_freq_offset import calc_freq_offset
from .calc_freq_offset import calc_freq_offset
#from ..tools.cassini_blocked import cassini_blocked
from ..tools.search_for_file import search_for_file
from ..tools.write_output_files import construct_filepath

import sys
sys.path.append('../../')
import rss_ringoccs as rss
sys.path.remove('../../')


class FreqOffsetFit(object):
    """Class to make a fit to extracted frequency offset. Uses predicted
    sky frequency, reconstructed sky frequency, and a fit to residual sky
    frequency to do so.
    """

    def __init__(self, rsr_inst, geo_inst, poly_order=7,
            f_uso_x=8427222034.34050, verbose=False):
        """
        Make a fit to sigma-clipped frequency offset.

        Args:
            rsr_inst
            geo_inst
        """

        # Check inputs for validity
        if not isinstance(rsr_inst, rss.rsr_reader.RSRReader):
            sys.exit('ERROR (FreqOffsetFit): rsr_inst input must be an '
                + 'instance of the RSRReader class')

        if not isinstance(geo_inst, rss.occgeo.Geometry):
            sys.exit('ERROR (FreqOffsetFIt): geo_inst input must be an '
                + 'instance of the Geometry class')

        if not isinstance(poly_order, int):
            print('WARNING (FreqOffsetFit): poly_order input must be an int. '
                + 'Ignoring current input and setting to order 9')
            poly_order = 9

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
            print('WARNING (freq_offset_fit.py): Invalid frequency band!')
            sys.exit()

        # Compute spline coefficients relating SPM to rho
        rho_geo_spl_coef = splrep(spm_geo, rho_geo)

        ### compute max and min SPM values for occultation ###
        # Evaluate spm-to-rho spline at raw SPM to get raw rho sampling
        #    that matches SPM values
        self.raw_rho = splev(self.raw_spm_vals,rho_geo_spl_coef)
        # Create boolean mask where True is within occultation range and
        #    False is outside the occultation range -- this is generalized
        #    to work for diametric and chord occultations
        inds = [(self.raw_rho>6.25e4)&(self.raw_rho<3e5)]
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
        if verbose:
            print('\tCalculating predicted frequency offset...')
        spm0, f_sky_pred = rsr_inst.get_f_sky_pred(f_spm=f_spm)
        f_sky_recon = calc_f_sky_recon(f_spm, rsr_inst, sc_name, f_uso,
                kernels)

        # Interpolate rho to frequency time values
        f_rho = splev(f_spm, rho_geo_spl_coef)

        # Compute residual sky frequency
        f_sky_resid = f_offset - (f_sky_recon - f_sky_pred)
        if verbose:
            print('\tCreating sigma clipping mask array...')
        self.__fsr_mask = self.create_mask(f_spm, f_rho, f_sky_resid)

        # Fit frequency offset residual
        if verbose:
            print('\tCalculating fit to frequency offset residuals...')
        f_sky_resid_fit = self.fit_f_sky_resid(f_spm, f_rho, f_sky_resid,
                poly_order=poly_order)

        # Draw and save reference plot
        self.plotFORFit(f_spm,f_sky_resid,f_sky_resid_fit,self.__fsr_mask,
                        spm_min,spm_max,poly_order)

        # Calculate frequency offset fit
        self.f_offset_fit = f_sky_resid_fit + (f_sky_recon - f_sky_pred)
        self.f_spm = f_spm
        self.f_sky_pred  = f_sky_pred
        self.f_sky_resid_fit = f_sky_resid_fit

    def create_mask(self, f_spm, f_rho, f_sky_resid):
        """
        Purpose:
            Creates a Boolean mask array which excludes data based on the
            following critera:
                #. ring or planetary occultation in region prevents accurate
                   estimation of the offset frequency
                #. offset frequencies fall more than 5-sigma beyond the median
                   offset frequency
                #. adjacent data all excluded by previous requirements (excludes
                   noise which by happenstance satisfies the above criteria)

        Arguments:
            :f_spm (*np.ndarray*): SPM sampled by ``calc_freq_offset`` when calculating
                                    the offset frequencies for the occultation
            :f_rho (*np.ndarray*): ring intercept radius of the spacecraft signal
                                    resampled to match f_spm
            :f_sky_resid (*np.ndarray*): residual sky frequency
        """

        # Create mask array that includes everything
        fsr_mask = np.array([True for i in range(len(f_sky_resid))],dtype=bool)

        # Compute median, standard deviation, and implememt sigma-clipping
        #   for data which fall in acceptable regions
        fsr_median = np.nanmedian(f_sky_resid[fsr_mask])
        fsr_stdev = 3.*np.sqrt(np.nanmedian(np.square(f_sky_resid-fsr_median)))
        if fsr_stdev < 1 or fsr_stdev > 20 :
            fsr_stdev = 1

        # iteratively check to see if each residual value is within 3 sigma
        for i in range(len(f_sky_resid)):
            if (f_sky_resid[i] < fsr_median - fsr_stdev) or (f_sky_resid[i] > fsr_median + fsr_stdev):
                fsr_mask[i] = False

        ## Polynomial fit clipping
        # try a 9th order polynomial fit
        pinit = np.polyfit(f_spm[fsr_mask], f_sky_resid[fsr_mask], 9)#np.polyfit(f_spm[fsr_mask], f_sky_resid[fsr_mask], 9)
        # Compute standard deviation from fit and implememt sigma-clipping
        #   for data which fall in acceptable regions
        fit_stdev = 3.*np.sqrt(np.nanmedian(np.square(f_sky_resid-np.polyval(pinit,f_spm))))
        # if the fit can give us a reasonable constraint, use it to help sigma clip
        if fit_stdev < 2 :
            # Create new mask array that includes everything
            fsr_mask = np.array([True for i in range(len(f_sky_resid))],dtype=bool)
            # iteratively check to see if each residual value is within 3 sigma
            for i in range(len(f_sky_resid)):
                if (f_sky_resid[i] < np.polyval(pinit,f_spm[i]) - fit_stdev) or (f_sky_resid[i] > np.polyval(pinit,f_spm[i]) + fit_stdev):
                    fsr_mask[i] = False
        else:
            fit_stdev = 0.1

        ## iteratively check adjacent values for false positives -- i.e.,
        #       all four adjacent mask array values are False
        #       first forwards
        for i in range(2,len(fsr_mask)-2):
            if fsr_mask[i]:
                if not fsr_mask[i-2]:
                    if not fsr_mask[i-1]:
                        if not fsr_mask[i+1]:
                            if not fsr_mask[1+2]:
                                fsr_mask[i] = False
        # now check backwards, just in case false positives were supporting
        # each other and preventing removal
        for i in range(len(fsr_mask)-2,2,-1):
            if fsr_mask[i]:
                if not fsr_mask[i-2]:
                    if not fsr_mask[i-1]:
                        if not fsr_mask[i+1]:
                            if not fsr_mask[1+2]:
                                fsr_mask[i] = False

        ## return frequency sky residual mask array
        return fsr_mask

    def fit_f_sky_resid(self, f_spm, f_rho, f_sky_resid,poly_order=None, verbose=False):
        """
        Fit a polynomial to residual frequency.

        Arguments:
            :f_spm (*np.ndarray*): SPM sampled by ``calc_freq_offset`` when calculating
                                    the offset frequencies for the occultation
            :f_rho (*np.ndarray*): ring intercept radius of the spacecraft signal
                                    resampled to match f_spm
            :f_sky_resid (*np.ndarray*): residual sky frequency

        Keyword Arguments:
            :poly_order (*float*): Order of polynomial fit to residual frequency
            :verbose (*bool*): If True, print processing steps
        """

        if not isinstance(poly_order, int):
            print('WARNING (FreqOffsetFit): poly_order input must be an int. '
                + 'Ignoring current input and setting to order 9')
            poly_order = 9


        npts = len(f_spm)
        spm_temp = ((f_spm - f_spm[int(npts / 2)])
            / max(f_spm - f_spm[int(npts / 2)]))

        ## fit using polynomial of user-selected order
        coef = np.polyfit(spm_temp[self.__fsr_mask],f_sky_resid[self.__fsr_mask],
                                poly_order)

        '''if verbose:
            print('\tPolynomial sum squared residuals:',stats[0])'''

        f_sky_resid_fit = np.polyval( coef, spm_temp )

        return f_sky_resid_fit

    # Create and save a plot of the offset residual fit
    def plotFORFit(self,spm,resid,fit,mask,spm_min,spm_max,poly_order):
        """
        Purpose:
            Plot results of the automated frequency offset residual fit and save
            plot to a file. File name will match the *.LBL and *.TAB nomenclature.

        Arguments:
            :spm (*np.ndarray*): SPM sampled by ``calc_freq_offset`` when calculating
                                    the offset frequencies for the occultation
            :resid (*np.ndarray*): residual sky frequency
            :fit (*np.ndarray*): polynomial fit to the residual sky frequency
            :mask (*np.ndarray*): boolean array used to mask residual sky frequency
                                    for the polynomial fitting
            :spm_min (*float*): start of occultation in SPM
            :spm_max (*float*): end of occultation in SPM
            :poly_order (*float*): order of polynomial fit to the residual sky frequency
        """
        # residuals used for fit
        plt.plot(spm[mask],resid[mask],'.k')
        # all residuals
        plt.plot(spm,resid,'-',color='0.5',lw=1)
        # indicate limits for ring system
        linds = [(self.raw_rho>=7e4)&(self.raw_rho<=1.4e5)]
        imin = np.argmin(self.raw_rho>=7.4e4)
        imax = np.argmax(self.raw_rho<=1.4e5)
        plt.axvline(self.raw_spm_vals[imin],dashes=[12,4],color='0.2')
        plt.axvline(self.raw_spm_vals[imax],dashes=[12,4],color='0.2')
        # fit to residuals
        plt.plot(spm,fit,'-r')
        # limits to plot
        plt.xlim(spm_min-100,spm_max+100)
        plt.ylim(np.nanmin(resid[mask])-0.1,np.nanmax(resid[mask])+0.1)
        # labels
        plt.xlabel('SPM (sec)')
        plt.ylabel(r'$f_{predict}-f_{observe}$')
        plt.title('Frequency Offset Residual Fit for PolyOrder '+str(poly_order))
        #generate plot file names
        filenames,outdirs = construct_filepath(self.rev_info,'FORFIT')
        # output
        for file,dir in zip(filenames,outdirs):
            plt.savefig(dir+file+'.PDF')
        plt.close()
"""
Revisions:
    2018 Oct 11 - jfong - copied from v1.0 freq_offset_fit.py
    2018 Nov 15 - sflury - changed sigma clipping so that minimum standard dev
                           is 1 Hz to prevent exclusion of actual trend
                         - added sigma clipping based on initial polynomial fit
                           in creat_mask to improve data selection
                         - tweaked plotting output for better visualization
"""
