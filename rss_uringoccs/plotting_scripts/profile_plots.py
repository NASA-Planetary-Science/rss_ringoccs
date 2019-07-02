### BEGIN USER INPUT
date = '20190522'   # processing date
sn = '0001'         # processing serial number
### END USER INPUT
import os
import numpy as np
import matplotlib.pyplot as plt
### suppress any annoying warning messages
import warnings as wrn # for error suppression
wrn.filterwarnings('ignore')
np.seterr(all='ignore')
# create directory structure if not in place
if not os.path.isdir('../output/PLOTS'):
    os.mkdir('../output/PLOTS')
if not os.path.isdir('../output/PLOTS/DLP005M'):
    os.mkdir('../output/PLOTS/DLP005M')
if not os.path.isdir('../output/PLOTS/TAU020M'):
    os.mkdir('../output/PLOTS/TAU020M')
if not os.path.isdir('../output/PLOTS/TAU050M'):
    os.mkdir('../output/PLOTS/TAU050M')
# ring name and location data
lab = ['6','5','4',r'$\alpha$',r'$\beta$',r'$\eta$',r'$\gamma$',
        r'$\delta$',r'$\varepsilon$']
# dictionary of ring locations in SPM
rings_km = {
    'E':{'6':41794.93, '5':42155.17, '4':42555.67,
         'A':44686.59, 'B':45673.33, 'N':47176.51,
         'G':47628.26, 'D':48297.35, 'E':51342.23},
    'I':{'6':41871.35, '5':42304.04, '4':42556.28,
         'A':44736.75, 'B':45640.74, 'N':47176.28,
         'G':47621.59, 'D':48300.57, 'E':50796.85}}
# ring widths
width = {
    'I':{'6':1.52, '5':2.75, '4':1.95,'A':10.59, 'B':7.03, 'N':1.54,
         'G':3.83, 'D':6.7, 'E':22.43},
    'E':{'6':1.72, '5':2.62, '4':2.67, 'A':4.22, 'B':11.19, 'N':1.53,
         'G':1.63, 'D':2.7, 'E':74.93}}
# plot reconstructed profiles
for res in ['050','020']:
    for dir in ['I','E']:
        for rl,ring in zip(lab,['6','5','4','A','B','N','G','D','E']):

            rmin = rings_km[dir][ring]-width[dir][ring]
            rmax = rings_km[dir][ring]+width[dir][ring]

            subdir = '../output/'+dir+'/'+ring+'/'
            file = 'VGR2_X43_'+dir+'_URING_'+ring+'_TAU_00'+res+'M_'+date+'_'+sn+'.TAB'
            tau = np.loadtxt(subdir+file,delimiter=',').T

            tmax = np.max(tau[6][(tau[0]>rmin)&(tau[0]<rmax)])

            plt.plot(tau[0],tau[6],'k',lw=1)
            plt.plot(tau[0],tau[8],'k',dashes=[12,6])
            if width[dir][ring] > 10 :
                plt.xlim(rmin,rmax)
            else:
                plt.xlim(rings_km[dir][ring]-10,rings_km[dir][ring]+10)
            if tmax < np.min(tau[8])-0.25 :
                plt.ylim(np.max(tau[8])+0.25,-0.5)
            else:
                plt.ylim(tmax+0.25,-0.5)
            plt.xlabel(r'RING INTERCEPT RADIUS (KM)')
            plt.ylabel(r'NORMAL OPTICAL DEPTH')
            plt.title(dir+' '+rl+' RECONSTRUCTION AT '+res+' M')
            plt.tight_layout()
            plt.savefig('../output/PLOTS/TAU'+res+'M/'+file[:-18]+'.PNG')
            plt.close()
# plot diffraction-limited profiles
for dir in ['I','E']:
    for rl,ring in zip(lab,['6','5','4','A','B','N','G','D','E']):

        subdir = '../output/'+dir+'/'+ring+'/'
        file = 'VGR2_X43_'+dir+'_URING_'+ring+'_DLP_005M_'+date+'_'+sn+'.TAB'
        rho,pow,phase = np.loadtxt(subdir+file,delimiter=',',usecols=(0,5,7)).T

        n = int(100)
        psmooth = np.convolve(pow,np.ones(n)/float(n),'same')
        r = [(rho>rmin)&(rho<rmax)]

        fig,ax = plt.subplots(2,1,sharex=True,figsize=(6,6),gridspec_kw={'hspace':0})
        ax[1].plot(rho,pow,color='0.5',lw=1,label='raw')
        ax[1].plot(rho,psmooth,color='0.0',lw=1,label='smoothed')
        ax[1].set_ylim(0,2)
        if width[dir][ring] > 10 :
            plt.xlim(rmin,rmax)
        else:
            plt.xlim(rings_km[dir][ring]-25,rings_km[dir][ring]+25)
        ax[1].set_xlabel(r'RING INTERCEPT RADIUS (KM)')
        ax[1].set_ylabel(r'NORMALIZED POWER')

        psmooth = np.convolve(phase,np.ones(n)/float(n),'same')
        ax[0].plot(rho,phase/np.pi,color='0.5',lw=1,label='raw')
        ax[0].plot(rho,psmooth/np.pi,color='0.0',lw=1,label='smoothed')
        ax[0].set_ylim(-0.5,0.5)
        ax[0].set_ylabel(r'PHASE (CYCLES SEC$^{-1}$)')
        ax[0].set_title(dir+' '+rl+' DIFFRACTION-LIMITED PROFILE AT 005 M')
        plt.tight_layout()
        plt.savefig('../output/PLOTS/DLP005M/'+file[:-4]+'.PNG')
        plt.close()
