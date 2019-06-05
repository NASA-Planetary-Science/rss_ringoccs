import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../')
from rss_ringoccs.scatter.spectro_reader import read_spectro
sys.path.remove('../')
import os

mpath = '../output/'

paths = [mpath+dir+'/'+subdir+'/'
            for dir in os.listdir(mpath)
            if '.' not in dir and float(dir[3:]) < 8
            for subdir in os.listdir(mpath+dir)
            if '.' not in subdir]

for path in paths:

    for dir in os.listdir(path):

        subdir = path + dir

        files = os.listdir(subdir)
        sfile = [file for file in files if '_SCATTER_' in file and file.endswith('.TAB')][0]
        tfile = [file for file in files if '_TAU_' in file and file.endswith('.TAB')][0]

        t,r,f,S = read_spectro(subdir+'/'+sfile)

        rho,tau = np.loadtxt(subdir+'/'+tfile,delimiter=',',usecols=(0,6)).T


        # frequency filtering to include only the thermal receiver power,
        #   averaged over time
        S_freq_filt = np.nanmean(S[((f>-450)&(f<-200))|((f>200)&(f<450)),:],0)
        # time filtering to include only the signal outside the occultation
        # by looking at just the first and last 1,000 seconds
        rmin = 7.3e4
        rmax = 1.39e5
        p_noise = np.nanmean(S_freq_filt[(r<rmin)|(r>rmax)])
        vmin = p_noise#/2

        for i in range(len(t)):
            S[(f>-15)&(f<15),i] = vmin
            S[(S[:,i]>0.3),i] = vmin

        print(subdir)
        print(sfile)




        fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,4),dpi=128)
        # plot profile
        ax[0].plot(rho/1e3,tau,'-k',lw=1)
        ax[0].tick_params(axis='x',direction='in',which='both',
                            top=True,bottom=True,labelbottom=False,
                            labeltop=False,length=8)
        ax[0].set_ylim(-0.1,4)
        #ax[0].set_ylabel('Normal Optical Depth')
        fig.text(0.05,0.875,'Normal Optical Depth',rotation=90)
        ax[0].set_title(path.split('/')[3]+' '+sfile[13:16])#[0:18])
        # plot spectrogram
        im = ax[1].pcolormesh(r/1e3,f/1e3,10.*np.log10(S)-10,
                                vmin=10.*np.log10(vmin)-10,cmap='plasma')
        # plot ticks, limits, and labels
        ax[1].tick_params(axis='x',direction='out',which='both',
                            top=False,bottom=True,labelbottom=True,
                            labeltop=False,length=6)
        ax[1].set_xlabel('$\\rho$ ($10^3$ km)')
        ax[1].set_ylim(-0.45,0.45)
        ax[1].set_ylabel('Frequency (kHz)')
        cbar_ax = fig.add_axes([0.9, 0.15, 0.01, 0.73])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label(r'$P/P_0$ (dB)')
        plt.gcf().subplots_adjust(bottom=0.15)

        plt.subplots_adjust(wspace=0.0,hspace=0.0)
        ### Save different slices
        for rlim in np.arange(70,140,5):
            ax[0].set_xlim(rlim,rlim+5)
            if rlim < 100 :
                rtag = '0' + str(int(rlim))
            else:
                rtag = str(int(rlim))
            plt.savefig(subdir+'/'+sfile[0:18]+'_SCATGAL_RHO'+rtag+'_'+sfile[-17:-4]+'.PNG')
        ax[0].set_xlim(117.5,118.5)
        plt.savefig(subdir+'/'+sfile[0:18]+'_SCATGAL_RHO117_'+sfile[-17:-4]+'.PDF')

        plt.close()
        # requires ImageMagick!!
        '''
        os.system('convert '+subdir+'/'+sfile[0:18]+'*_SCATGAL_RHO*'+sfile[-17:-4]+'.PNG '+subdir+'/'+sfile[0:18]+'_SCATGAL_'+sfile[-17:-4]+'.PDF')

        for rlim in np.arange(70,140,5):
            ax[0].set_xlim(rlim,rlim+5)
            if rlim < 100 :
                rtag = '0' + str(int(rlim))
            else:
                rtag = str(int(rlim))
            os.unlink(subdir+'/'+sfile[0:18]+'_SCATGAL_RHO'+rtag+'_'+sfile[-17:-4]+'.PNG')
        '''
