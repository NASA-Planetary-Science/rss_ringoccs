#!/usr/bin/env python
"""
Purpose:
    Normalize frequency-corrected power using a polynomial fit of specified order.

"""

import numpy as np
from scipy import signal
from scipy.interpolate import splrep
from scipy.interpolate import splev
import sys
sys.path.append('../..')
import rss_ringoccs as rss
sys.path.remove('../..')

from ..tools.cassini_blocked import cassini_blocked
from .namegen import plotname
import pdb
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


class Normalization(object):
    """
    Purpose:
        Finds freespace power based on ring geometry, locations of gaps
        computed by the Geometry class, and signal power relative to median
        power within the gaps. Fits the freespace power with a polynomial.
        If desired, this fitting process can be interactive. A plot will be
        saved following rss_ringoccs nomenclature to visualize the fit results.

    Arguments:
        :spm_raw (*np.ndarray*): SPM as sampled in the raw data
        :IQ_c (*np.ndarray*): frequency-corrected signal corresponding to spm_raw
        :geo_inst (*np.ndarray*):
                       instance of the Geometry class, used to estimate
                       the freespace regions within and surrounding the
                       ring system, accounting for Saturn occultation
                       and the Keplerian geometry of the rings
        :rho_km_vals (*np.ndarray*):
                       radial intercept poin of occultation in km

    Keyword Arguments
        :verbose (*bool*): when True, outputs information to
                       command line about the freespace power fitting
                       progress and results. All results will be output
                       to the *CAL*.LBL file regardless of this keyword.
                       Default is False.
        :order (*float*): a whole number specifying the order of the polynomial
                        fit to the freespace power. Default is 3.
        :interact (*bool*): If True, allows user to interactively adjust fit
                        to the freespace power. Default is False.
    """

    def __init__(self, spm_raw, IQ_c, geo_inst, rsr_inst, order=3, fittype='poly',
                                                interact=False, verbose=False):

        # store rev info
        self.rev_info = geo_inst.rev_info

        # Resample geometry values to raw sampling rate
        rho_km_vals = geo_inst.rho_km_vals
        #spm_to_rho = splrep(geo_inst.t_oet_spm_vals, geo_inst.rho_km_vals)
        #rho_km_full = splev(rsr_inst.spm_vals, spm_to_rho)

        #self.rho_to_spm = splrep(geo_inst.rho_km_vals, geo_inst.t_oet_spm_vals)
        freespace_spm = geo_inst.freespace_spm

        # downsample IQ so that diffraction fringes do not affect fit
        spm_down, p_obs_down = self.downsample_IQ(spm_raw, IQ_c,dt_down=0.2)
        self.spm_down = spm_down

        #rho_down = splev(spm_down, spm_to_rho)
        spm_to_rho = splrep(geo_inst.t_oet_spm_vals, geo_inst.rho_km_vals)
        rho_down = splev(spm_down, spm_to_rho)

        spm_to_rimp = splrep(geo_inst.t_oet_spm_vals, geo_inst.R_imp_km_vals)
        rimp_down = splev(spm_down, spm_to_rimp)

        # Create mask array based on freespace region predictions
        self.create_mask(spm_down, freespace_spm,p_obs_down)

        # Compute fit
        self.fit_freespace_power(spm_down, p_obs_down, order=order,fittype=fittype)
        # get User input to see if everything looks okay or needs tweaking
        #    will return modified freespace regions
        if interact :
            new_gaps = self.fit_check(spm_down,p_obs_down,freespace_spm,order)
        else:
            new_gaps = freespace_spm
        # plot final fit for user reference
        self.plot_power_profile(spm_down,p_obs_down,new_gaps,order,save=True)


    def fit_check(self,spm_down,p_obs_down,freespace_spm,order):
        """
        Purpose:
            Allows user to update the freespace regions and fit order during the
            freespace power fitting step. This is done by prompting the user for
            input in the command line and displaying the results of their input
            for the polynomial fit to the freespace power.
            Only called if the Normalization keyword ``interact`` is set to True.

        Arguments:
            :gaps_str (*str*): string containing the user input freespace regions.

        Returns:
            :gaps (*list*): an Nx2 list of floats indicating the lower and upper
                            limits to the user-specific freespace regions.
        """
        # plot fit for User to see
        self.plot_power_profile(spm_down,p_obs_down,freespace_spm,order)

        # Prompt if fit looks okay
        cont = 'y'#input('\nDo you want to continue with this fit? (y/n): ')

        new_gaps = freespace_spm
        while 'n' in cont or 'N' in cont:
            # Prompt if freespace regions to change
            change_freespace = input('\nDo you want to change the freespace regions? ')
            if 'y' in change_freespace or 'Y' in change_freespace:
                # print most recent freespace regions used
                print('\n')
                print('Last used freespace gaps in SPM: ')
                print('['+',\n'.join(str(g) for g in new_gaps)+']')
                print('\n')
                # Prompt for reversion to default (in case of oopses)
                revert_gaps = input('Do you want to revert to the default freespace gaps? ')
                # if reverting, then revert
                if 'y' in revert_gaps or 'Y' in revert_gaps:
                    new_gaps = freespace_spm
                # if not reverting, prompt for new freespace regions
                else:
                    new_gaps_str = input('Please input new freespace gaps in SPM and '
                        + 'press enter twice: ')
                    print('\n')
                    # if the string doesn't start with double brackets (not an NxM list)
                    #    then keep current gaps
                    if '[[' not in new_gaps_str and ']]' not in new_gaps_str :
                        print('Invalid entry for freespace gaps.\nMust be Nx2 list.')
                        print('Reverting to previous freespace regions.')
                        new_gaps = self.gaps
                    # otherwise, extract new freespace regions
                    else:
                        new_gaps = self.extract_list_from_str(new_gaps_str[1:-1])
                # create new mask based on changes to freespace regions
                self.create_mask(spm_down, new_gaps, p_obs_down/np.nanmax(p_obs_down))

            # Prompt for fit type
            print('\n')
            fittype = input('What fit type would you like to use? ' +
                    '(poly/spline): ')
            # Prompt for fit order
            print('\n')
            order = int(input('What fit order would you like to use?: '))
            print('\n')

            # update fit
            self.fit_freespace_power(spm_down, p_obs_down, order=order,fittype=fittype)
            # update plot for user
            self.plot_power_profile(spm_down,p_obs_down,new_gaps,order)
            cont = input('Do you want to continue with this fit? (y/n): ')

        return new_gaps


    def extract_list_from_str(self, gaps_str):
        """
        Purpose:
            Extract an Nx2 list from the string of user input freespace regions.

        Arguments:
            :gaps_str (*str*): string containing the user input freespace regions.

        Returns:
            :gaps (*list*): an Nx2 list of floats indicating the lower and upper
                            limits to the user-specific freespace regions.
        """
        gaps = []

        ind = [i for i, ltr in enumerate(gaps_str) if ltr=='[']

        for i in range(len(ind)):
            if i == len(ind)-1:
                str1 = gaps_str[ind[i]+1:-1]
            else:
                s1 = ind[i]+1
                s2 = ind[i+1]-3
                str1 = gaps_str[s1:s2]#ind[0]+1:ind[1]-3]
            f1 = float(str1.split(',')[0])
            f2 = float(str1.split(',')[1])
            print([f1,f2])
            gaps.append([f1,f2])



        return gaps

    def create_mask(self, spm_full, gaps_spm, pc):
        """
        Arguments:
            :geo_inst (*object*):
                           Instance of the Geometry class, used to estimate
                           the freespace regions within and surrounding the
                           ring system, accounting for Saturn occultation
                           and the Keplerian geometry of the rings
            :rho_km_vals (*np.ndarray*): radial intercept poin of occultation in km

        Sets attributes
            :__mask (*np.ndarray*):
                           array of booleans wherein True designates data
                           corresponding to freespace power and False
                           data corresponding to occultation events
            :__gaps (*list*):
                           an Nx2 list of lower and upper radial limits in
                           km for each of N gaps designated
        """

        # normalize corrected power to max corrected power inside the occultation
        pc_max = np.nanmax(pc[(spm_full>=gaps_spm[1][1])&(spm_full<=gaps_spm[-2][0])])
        pc_norm = pc/pc_max

        # get lower, upper radial limits to planet/atmosphere occultation
        # get lower, upper radii for each gap in C-ring, Cassini Division,
        #    and the Enke gap

        # storage
        pc_median = []
        ### compute median power in each freespace region
        #for rho_limits in gaps:
        #    # Boolean mask including only rho values within gap
        #    ind = [(rho_km_full>rho_limits[0])&(rho_km_full<rho_limits[1])]
        #    # find median corrected power within space
        #    pc_median += [np.nanmedian(pc_norm[ind])]
        gaps_spm_copy = gaps_spm.copy()

        for spm_limits in gaps_spm_copy:
            # Boolean mask including only spm values within gap
            ind = [(spm_full>=spm_limits[0])&(spm_full<=spm_limits[1])]
            # find median corrected power within space
            #if pc_median < np.median(pc_median)*1.3 and pc_median > np.median(pc_median)*1.3:
            pcm = np.nanmedian(pc_norm[ind])
            if pcm < 0.5 or pcm > 1.25 :
                gaps_spm.remove(spm_limits)
            else:
                pc_median += [pcm]#[np.nanmedian(pc_norm[ind])]
            #else:
            #    gaps_spm.remove(spm_limits)

        # create mask array that includes only freespace power
        fsp_mask = np.array([False for i in range(len(spm_full))])

        # iterate over all radii
        for i in range(len(spm_full)):
            # check to see if intercept radius falls within a gap
            #for rho_limits in gaps:
            for spm_limits,p_c_m in zip(gaps_spm,pc_median):
                # looks at specific gap and compares radius to gap boundaries
                if spm_full[i] >= spm_limits[0] and spm_full[i] <= spm_limits[1]:
                    #fsp_mask[i] = True
                    # compares the corrected, normalized IQ to the median free-
                    # space power within the gap, TRUE if within 10%
                    if pc_norm[i] > p_c_m-0.1 and pc_norm[i] < p_c_m+0.1 :
                        # if near median power, and within gap, then change mask
                        # from exclude to include
                        fsp_mask[i] = True


        self.mask = fsp_mask
        self.gaps = gaps_spm


    def fit_freespace_power(self,spm,power,order=3,fittype='poly'):
        """
        Arguments:
            :spm (*np.ndarray*): downsampled SPM
            :power (*np.ndarray*): absolute square of downsapmled phase-corrected signal

        Keyword Arguments:
            :order (*float*): order of the fit, whole number between 1 and 5.
                                Default order is 3.
            :type (*str*): type of fit to use, default is 'poly'. Options are
                            - 'poly' a single polynomial
                            - 'spline' an unsmoothed spline fit

        Returns
            :fit (*np.ndarray*): best fit to freespace power
        """
        self.fittype=fittype
        self.order=order

        # polynomial fit
        if fittype == 'poly' :

            # compute polynomial coefficients for a least squares fit to the
            #   power in the given freespace regions
            coef = np.polyfit(spm[self.mask],power[self.mask],order)
            # evaluate the polynomial at all spm
            fit = np.polyval(coef, spm)
            chi2 = np.sum(np.square(np.polyval(coef,spm[self.mask])-power[self.mask]))

            #fit = np.polyval(coef, spm_full)

        # spline fit
        elif fittype == 'spline' :

            # compute knots as midpoint between gap boundaries
            knots_spm = []
            for i in range(len(self.gaps)):
                #knots+=[np.nanmean(self.__gaps[i])]
                knots_spm+=[np.nanmean(self.gaps[i])]

            #knots_spm = splev(knots, self.rho_to_spm)
            mask = self.mask

            ml = len(self.mask)
            smooth = (ml-np.sqrt(2*ml),ml+np.sqrt(2*ml))

            # use knots and given order to compute B-spline coefficients
            coef = splrep(spm[self.mask],power[self.mask],k=order,t=knots_spm, s=smooth)

            # evaluate spline at all spm
            fit = splev(spm,coef)
            chi2 = np.sum(np.square(splev(spm[self.mask],coef)-power[self.mask]))

        else:
            # if keyword type not recognized, state for user and fit with poly
            print('Unknown option ',fittype,' given for keyword \'type\'')
            print('Proceeding with a polynomial fit')
            # compute polynomial coefficients for a least squares fit to the
            #   power in the given freespace regions
            coef = np.polyfit(spm[self.mask],power[self.mask],order)
            # evaluate the polynomial at all spm
            fit = np.polyval(coef, spm)
            chi2 = np.sum(np.square(np.polyval(coef,spm[self.mask])-power[self.mask]))


        self.pnorm_fit = fit
        self.chi_squared = chi2
        return None

    def plot_power_profile(self,spm,pow,gaps,order,save=False):
        """
        Purpose:
            Plot results of the freespace power for total profile and each individual
            freespace region and either show plot in a GUI or save it to a file.
            File name will match the *.LBL and *.TAB nomenclature.

        Arguments:
            :spm (*np.ndarray*): SPM sampled by ``calc_freq_offset`` when calculating
                                    the offset frequencies for the occultation
            :pow (*np.ndarray*): residual sky frequency
            :gaps (*np.ndarray*): gap edges used to select freespace power
            :order (*float*): order of polynomial fit to the residual sky frequency

        Keyword Arguments:
            :save (*bool*): If True, saves plots to file. Otherwise, plot is shown in GUI.
                            Default is False.
        """
        # max power inside the occultation
        pmax = np.nanmax(pow[(spm>=gaps[0][0])&(spm<=gaps[-1][1])])
        # scaling of SPM for visualization
        spm_off = 10. * round(spm[0]/1e4)

        # find number of columns to make gap plots
        rmdr = float(len(gaps))%4
        if rmdr > 0 :
            ncols = int(len(gaps)/4)+1
        else :
            ncols = int(len(gaps)/4)

        # draw figure and create subplots
        plt.close('all')
        fig = plt.figure(figsize=(1+1.25*float(ncols),8))
        gs = GridSpec(6,ncols)
        ### total profile - takes two rows ###
        ax_prof = fig.add_subplot(gs[0:2,:])
        ax_prof.plot(spm[self.mask]/1e3-spm_off, pow[self.mask]/pmax, '.k',ms=6,zorder=1)
        ax_prof.plot(spm/1e3-spm_off, pow/pmax, color='blue',zorder=2)
        ax_prof.plot(spm/1e3-spm_off, self.pnorm_fit/pmax, color='red',zorder=3)
        ax_prof.set_xlim((gaps[0][0]-100)/1e3-spm_off,(gaps[-1][1]+100)/1e3-spm_off)
        ax_prof.set_ylim(0,1.1)
        ### individual freespace regions ###
        # add subplots
        ax = []
        for i in range(4):
            for j in range(ncols):
                ax+=[fig.add_subplot(gs[i+2,j])]
        # plot regions
        for i in range(len(ax)):
            # check to make sure there's a gap to plot
            if i < len(gaps):
                # x range and limits
                dx = gaps[i][1]-gaps[i][0]
                xmin = gaps[i][0]-5#0.5*dx
                xmax = gaps[i][1]+5#0.5*dx
                # label on profile
                lims = np.array(gaps[i])/1e3-spm_off
                c=np.random.rand(3)
                #ax_prof.text(np.mean(lims),1.0,str(i+1),zorder=4)
                #ax_prof.fill_between(lims,-1,2,color=c,zorder=0)
                # plotting
                #ax[i].fill_between(lims,-1,2,color=c,zorder=0)
                ax[i].plot(spm[self.mask]/1e3-spm_off, pow[self.mask]/pmax,'.k',ms=6)
                ax[i].plot(spm/1e3-spm_off, pow/pmax, color='blue')
                ax[i].plot(spm/1e3-spm_off, self.pnorm_fit/pmax, color='red')

                ax[i].set_xlim(xmin/1e3-spm_off,xmax/1e3-spm_off)
                ax[i].set_ylim(0,1.1)
                # remove tick labels for land-locked plots
                if i%ncols > 0 :
                    for tck in ax[i].get_yticklabels():
                        tck.set_visible(False)
            # if not, then clear the plot frame
            else:
                ax[i].axis('off')

        # subplot margin setting
        plt.subplots_adjust(wspace=None, hspace=0.3,left=0.1,right=0.975,top=0.975,bottom=0.075)
        # axis labels
        fig.text(0.45,0.01,r'SPM - '+str(int(spm_off))+' ($10^3$ sec)')
        fig.text(0.01,0.5,r'Power (arb.)',rotation=90)
        if save:
            filename = plotname(self.rev_info,'FSPFIT')
            plt.savefig(filename,dpi=128)
            plt.close('all')
        else:
            plt.show(block=False)

    def downsample_IQ(self, spm_raw, IQ_c_raw, dt_down=0.5, verbose=False):
        """
        Purpose:
            Downsample complex signal to specified time spacing to avoid
            diffraction pattern affecting freespace power fit

        Arguments:
            :spm_raw (*np.ndarray*): raw SPM values
            :IQ_c_raw (*np.nparray*): :math:`I_c+iQ_c` sampled at the raw SPM rate

        Keyword Arguments:
            :dt_down (*float*):
                Time spacing to downsample to
            :verbose (*bool*):
                If True, prints downsampled results

        Returns:
            :spm_vals_down (*np.ndarray*):
                SPM values after downsampling
            :rho_km_vals_down (*np.ndarray*):
                Rho values after downsampling
            :p_obs_down (*np.ndarray*):
                Observed power after downsampling
        """

        # Downsampling coefficient q
        dt_raw = spm_raw[1] - spm_raw[0]
        q = round(dt_down / dt_raw)

        # Downsample IQ_c by factor of q and not power because this is to
        #     match the resampling done in norm_diff_class.py, where it also
        #     resamples IQ_c
        IQ_c_down = signal.resample_poly(IQ_c_raw, 1, q)
        p_obs_down = abs(IQ_c_down)**2

        # New SPM, rho, and p_obs, at the downsampled resolution
        spm_vals_down = np.linspace(spm_raw[0], spm_raw[-1],
            num=len(p_obs_down))
        #rho_km_vals_down = self.__rho_interp_func(spm_vals_down)

        return spm_vals_down, p_obs_down#,rho_km_vals_down, p_obs_down
"""
Revisions:
      power_normalization.py
      2018 Oct 17 - jfong - edit inputs to create_mask
                            - remove rho_km_vals because already in geo_inst
                            - remove spm_vals -- never used?
                          - use atmos_occ_spm_vals in place of cassini_blocked
                          - index ind not m for IQ_c_median
                          - IQ_c_median is at raw spm values, not same dim
      2018 Nov 12 - sflury - changed median freespace power criterion from 1e-3
                             to 0.5 in the create_mask function
                           - moved power normalization plot to own function
                           - moved iterative user interaction to own function
      2018 Nov 15 - sflury - added clipping of freespace power that is too high,
                             made interactive mode optional with kwarg "interact"
"""
