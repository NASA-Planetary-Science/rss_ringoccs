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
if not os.path.isdir('../output/PLOTS/GALLERY'):
    os.mkdir('../output/PLOTS/GALLERY')
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
# plot all resolutions
for res in ['050','020']:

    for dir in ['I','E']:

        fig,axes = plt.subplots(3,3,figsize=(9,8))
        axes = axes.flatten()#[axes[i][j] for i in range(axes) for j in range(axes[i])]

        for rl,ring,ax in zip(lab,['6','5','4','A','B','N','G','D','E'],axes):

            rmin = rings_km[dir][ring]-width[dir][ring]
            rmax = rings_km[dir][ring]+width[dir][ring]

            subdir = '../output/'+dir+'/'+ring+'/'
            file = 'VGR2_X43_'+dir+'_URING_'+ring+'_TAU_00'+res+'M_'+date+'_'+sn+'.TAB'
            tau = np.loadtxt(subdir+file,delimiter=',').T

            tmax = np.max(tau[6][(tau[0]>rmin)&(tau[0]<rmax)])
            clp = [(tau[0]>rings_km[dir][ring]-0.75*width[dir][ring])&(tau[0]<rings_km[dir][ring]+0.75*width[dir][ring])]

            ax.plot(tau[0],tau[6],'k',lw=1)
            ax.plot(tau[0][clp],tau[8][clp],'k',lw=1)
            if width[dir][ring] > 10 :
                ax.set_xlim(rmin,rmax)
            else:
                ax.set_xlim(rings_km[dir][ring]-10,rings_km[dir][ring]+10)
            if tmax < np.min(tau[8])-0.25 :
                ax.set_ylim(np.max(tau[8])+0.25,-0.5)
            else:
                ax.set_ylim(tmax+0.25,-0.5)
            ax.text(0.1,0.3,rl,transform=ax.transAxes,fontsize=14)
        # axis labels
        axes[7].set_xlabel('RING INTERCEPT RADIUS (KM)')
        axes[3].set_ylabel('NORMAL OPTICAL DEPTH')
        if dir == 'I':
            direct = 'INGRESS'
        elif dir == 'E':
            direct = 'EGRESS'
        fig.suptitle('URANUS RINGS DURING '+direct+' RECONSTRUCTED AT '+res+' M')
        fig.tight_layout()
        fig.subplots_adjust(top=0.95)
        plt.savefig('../output/PLOTS/GALLERY/GALLERY_'+dir+'_'+res+'M.PDF',dpi=128)
        plt.close()
