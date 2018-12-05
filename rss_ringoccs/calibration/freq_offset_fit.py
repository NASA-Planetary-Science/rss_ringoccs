"""

freq_offset_fit_jol.py

Purpose: Makes a fit to the frequency offset made from freq_offset.py, using
         the frequency offset, predicted sky frequency, reconstructed sky
         frequency, and a fit to residual frequency.

WHERE TO GET NECESSARY INPUT:
    rsr_inst: Use an instance of the RSRReader class, found inside of
        rss_ringoccs/rsr_reader/rsr_reader.py
    geo_inst: Use an instance of the Geometry class, found inside of
        rss_ringoccs/occgeo/calc_occ_geometry.py

Revisions:
    2018 Oct 11 - jfong - copied from v1.0 freq_offset_fit.py
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
from .plotcal import plotCal
#from ..tools.cassini_blocked import cassini_blocked
from ..tools.search_for_file import search_for_file
from ..tools.write_intermediate_files import write_intermediate_files

import sys
sys.path.append('../../')
import rss_ringoccs as rss
sys.path.remove('../../')


class FreqOffsetFit(object):
    """Class to make a fit to extracted frequency offset. Uses predicted
    sky frequency, reconstructed sky frequency, and a fit to residual sky
    frequency to do so.
    """

    def __init__(self, rsr_inst, geo_inst, poly_order=9,
            f_uso_x=8427222034.34050, file_search=False, verbose=False):
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
        sc_name = geo_inst.history['Input Variables']['spacecraft']

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

        '''# Search for frequency offset fit (FOF) text file
        fof_file = search_for_file(self.year, self.doy,
                self.band, self.dsn, self.profdir, 'FOF')

        # If no file found, calculate frequency offset and save FOF file
        if (fof_file == 'N/A'):
            print('\tNo frequency offset file found...\n'
                    + '\tCalculating frequency offset...')
            f_spm, f_offset, f_offset_history = calc_freq_offset(
                    rsr_inst)
            write_intermediate_files(rsr_inst.year,
                    rsr_inst.doy, rsr_inst.band, rsr_inst.dsn, self.profdir,
                    'FOF', {'f_spm':f_spm, 'f_offset':f_offset})
        else:
            # Read in parameters from offset file
            print('\tExtracting frequency offset from:\n\t\t'
                    + '/'.join(fof_file.split('/')[0:5]) + '/\n\t\t\t'
                    + fof_file.split('/')[-1])
            freq_offset_file_vals = np.loadtxt(fof_file)
            f_spm = freq_offset_file_vals[:, 0]
            f_offset = freq_offset_file_vals[:, 1]'''

        print('\tCalculating observed frequency offset...')
        foff_inst = calc_freq_offset(rsr_inst)
        f_spm, f_offset = foff_inst.f_spm, foff_inst.f_offset

        print('\tCalculating predicted frequency offset...')
        spm0, f_sky_pred = rsr_inst.get_f_sky_pred(f_spm=f_spm)
        f_sky_recon = calc_f_sky_recon(f_spm, rsr_inst, sc_name, f_uso,
                kernels)

        # Interpolate rho to frequency time values
        rho_geo_spl_coef = splrep(spm_geo, rho_geo)
        f_rho = splev(f_spm, rho_geo_spl_coef)

        # Compute residual sky frequency
        f_sky_resid = f_offset - (f_sky_recon - f_sky_pred)

        print('\tCreating sigma clipping mask array...')
        self.__fsr_mask = self.__create_mask(f_rho, f_sky_resid)

        # Search for FRFP file
        # Fit frequency offset
        print('\tCalculating fit to frequency offset residuals...')
        f_sky_resid_fit = self.fit_f_sky_resid(f_spm, f_rho, f_sky_resid,
                poly_order=poly_order)

        # Calculate frequency offset fit
        self.f_offset_fit = f_sky_resid_fit + (f_sky_recon - f_sky_pred)
        self.f_spm = f_spm
        self.f_sky_pred  = f_sky_pred
        self.f_sky_resid_fit = f_sky_resid_fit

    def __create_mask(self, f_rho, f_sky_resid, polyclip=False):
        """
        Creates a Boolean mask array which excludes data based on the
        following critera:
            1) ring or planetary occultation in region prevents accurate
               estimation of the offset frequency
            2) offset frequencies fall more than 5-sigma beyond the median
               offset frequency
            3) adjacent data all excluded by previous requirements (excludes
               noise which by happenstance satisfies the above criteria)

        Args:
            polyclip (bool):
                Boolean indicating whether to perform sigma clipping based
                on an initial polynomial fit instead of only using the
                data median.
        """
        # TODO polyclip kwarg

        # Create mask array that includes everything
        fsr_mask = np.array([True for i in range(len(f_sky_resid))],
                dtype=bool)

        # Exclude regions within B ring with high optical depths and
        #   regions far outside of the rings, which prevent any useful
        #   assessment of the frequency offset
        rho_exclude = [[0, 70000], [91900, 94000], [98000, 118000],
                [194400, np.inf]]

        for i in range(len(f_rho)):
            for rl in rho_exclude:
                if f_rho[i] > rl[0] and f_rho[i] < rl[1]:
                    fsr_mask[i] = False

        # Compute median, standard deviation, and implememt sigma-clipping
        #   for data which fall in acceptable regions
        fsr_median = np.median(f_sky_resid[fsr_mask])
        fsr_stdev = 5. * np.sqrt(np.median(np.square(f_sky_resid[fsr_mask]
            - fsr_median)))

        for i in range(len(f_sky_resid)):
            if (f_sky_resid[i] < fsr_median - fsr_stdev) or (f_sky_resid[i] > fsr_median + fsr_stdev):
                fsr_mask[i] = False

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

        for i in range(len(fsr_mask)-2,2,-1):
            if fsr_mask[i]:
                if not fsr_mask[i-2]:
                    if not fsr_mask[i-1]:
                        if not fsr_mask[i+1]:
                            if not fsr_mask[1+2]:
                                fsr_mask[i] = False

        ## return frequency sky residual mask array
        return fsr_mask

    def fit_f_sky_resid(self, f_spm, f_rho, f_sky_resid,
            poly_order=None, verbose=False):
        """
        Fit a polynomial to residual frequency.

        Args:
            poly_order (int):
                Order of polynomial fit made to residual frequency
            verbose (bool):
                Print processing steps
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

        frfp = {'coef': coef}

        write_intermediate_files(self.year, self.doy, self.band, self.dsn,
                self.profdir, 'FRFP', frfp)

        '''if verbose:
            print('\tPolynomial sum squared residuals:',stats[0])'''

        f_sky_resid_fit = np.polyval( coef, spm_temp )
        #print(coef)

        plotCal(f_spm,f_sky_resid,f_sky_resid_fit,self.__fsr_mask,'FORFIT',
                self.rev_info,'SPM',r'$f_{predict}-f_{observe}$',
                'Frequency Offset Residual Fit for PolyOrder '+str(poly_order))

        return f_sky_resid_fit
