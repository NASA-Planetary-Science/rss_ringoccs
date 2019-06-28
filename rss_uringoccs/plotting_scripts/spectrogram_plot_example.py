### BEGIN USER INPUT
sdate = '20190522'   # scattered signal date
ssn = '0001'         # scattered signal serial number
tdate = '20190522'   # processing date
tsn = '0001'         # processing serial number
res = '050'         # reconstruction resolution
poly_order = 9      # polynomial order
### END USER INPUT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../pipeline/')
import pipeline_params as pipe
sys.path.remove('../pipeline/')
sys.path.append('../')
from rss_ringoccs.scatter.spectro_reader import read_spectro
sys.path.remove('../')
### suppress any annoying warning messages
import warnings as wrn # for error suppression
wrn.filterwarnings('ignore')
np.seterr(all='ignore')
# create directory structure if not in place
if not os.path.isdir('../output/PLOTS'):
    os.mkdir('../output/PLOTS')
# create directory structure if not in place
if not os.path.isdir('../output/PLOTS/SCATTER'):
    os.mkdir('../output/PLOTS/SCATTER')
# ring name and location data
lab = ['6','5','4',r'$\alpha$',r'$\beta$',r'$\eta$',r'$\gamma$',
        r'$\delta$',r'$\varepsilon$']

for dir in ['I','E']:
    for i,name in enumerate(pipe.rnames):
        # load in
        refdir = '../output/'+dir+'/'+name+'/'
        sfile = 'VGR2_X43_'+dir+'_URING_'+name+'_SCATTER_'+sdate+'_'+ssn+'.TAB'
        tfile = 'VGR2_X43_'+dir+'_URING_'+name+'_TAU_00'+res+'M_'+tdate+'_'+tsn+'.TAB'

        t,r,f,S = read_spectro(refdir+sfile)
        rho,tau = np.loadtxt(refdir+tfile,delimiter=',',usecols=(0,6)).T

        # frequency filtering to include only the thermal receiver power,
        #   averaged over time
        S_freq_filt = np.nanmean(S[((f>-450)&(f<-200))|((f>200)&(f<450)),:],0)
        # time filtering to include only the signal outside the occultation
        # by looking at just the edges of the occultation
        rmin = np.nanmedian(rho)/1e3-3
        rmax = np.nanmedian(rho)/1e3+3
        p_noise = np.nanmean(S_freq_filt[(r/1e3<rmin)|(r/1e3>rmax)])
        vmin = p_noise

        for ti in range(len(t)):
            S[(f>-15)&(f<15),ti] = vmin
            S[(S[:,ti]>0.075),ti] = vmin

        fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,4),dpi=128)
        # plot profile
        ax[0].plot(rho/1e3,tau,'-k',lw=1)
        ax[0].tick_params(axis='x',direction='in',which='both',
                            top=True,bottom=True,labelbottom=False,
                            labeltop=False,length=8)
        ax[0].set_ylim(-0.1,4)
        #ax[0].set_ylabel('Normal Optical Depth')
        fig.text(0.05,0.875,'Normal Optical Depth',rotation=90)
        ax[0].set_title('URING'+' '+lab[i]+' '+dir)
        # plot spectrogram
        im = ax[1].pcolormesh(r/1e3,f/1e3,10.*np.log10(S)-10,vmax=-25,
                                vmin=-30,cmap='plasma')#10.*np.log10(vmin)-10
        # show wavelength dependence
        m = -0.55
        r0 = np.median(rho)
        ax[1].axvline(51.30476,color='c',dashes=[4,12],lw=1)
        ax[1].axvline(51.37969,color='c',dashes=[4,12],lw=1)
        #ax[1].plot(r/1e3,m*(r-r0)/1e3-0.1,'-c',dashes=[6,12])
        #ax[1].plot(r/1e3,m*(r-r0)/1e3+0.1,'-c',dashes=[6,12])
        ax[1].tick_params(axis='x',direction='out',which='both',
                            top=False,bottom=True,labelbottom=True,
                            labeltop=False,length=6)
        ax[1].set_xlabel('$\\rho$ ($10^3$ km)')
        ax[1].set_ylim(-1.1,1.1)
        ax[1].set_xlim(np.median(r/1e3)-2,np.median(r/1e3)+2)
        ax[1].set_ylabel('Frequency (kHz)')
        cbar_ax = fig.add_axes([0.9, 0.15, 0.01, 0.73])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label(r'$P/P_0$ (dB)')
        plt.gcf().subplots_adjust(bottom=0.15)

        plt.subplots_adjust(wspace=0.0,hspace=0.0)

        plt.savefig('../output/PLOTS/SCATTER/'+sfile[:-14]+'_TAU.PNG')
        plt.close()

        fig,ax = plt.subplots(1,1,sharex=True,figsize=(8,8),dpi=128)
        # plot spectrogram
        im = ax.pcolormesh(r/1e3,f/1e3,10.*np.log10(S)-10,vmax=-25,
                                vmin=-30,cmap='plasma')#10.*np.log10(vmin)-10
        # show wavelength dependence
        if name == 'E' and dir == 'E' :
            m = -0.55
            r0 = 51304.76
            ax.plot(r/1e3,m*(r-r0)/1e3,'-c',dashes=[6,12])
            r0 = 51379.69
            ax.plot(r/1e3,m*(r-r0)/1e3,'-c',dashes=[6,12])
        ax.tick_params(axis='x',direction='out',which='both',
                            top=False,bottom=True,labelbottom=True,
                            labeltop=False,length=6)
        ax.set_xlabel('$\\rho$ ($10^3$ km)')
        ax.set_ylim(-0.8,0.8)
        ax.set_xlim(np.median(r/1e3)-2,np.median(r/1e3)+2)
        ax.set_ylabel('Frequency (kHz)')
        cbar_ax = fig.add_axes([0.9, 0.15, 0.01, 0.73])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label(r'$P/P_0$ (dB)')
        ax.set_title('URING'+' '+lab[i]+' '+dir)

        plt.gcf().subplots_adjust(bottom=0.15)

        plt.subplots_adjust(wspace=0.0,hspace=0.0)

        plt.savefig('../output/PLOTS/SCATTER/'+sfile[:-14]+'.PNG')
        plt.close()
