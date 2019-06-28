### BEGIN USER INPUT
res = '00050'       # processing resolution
date = '20190522'   # processing date
sn = '0002'         # processing serial number
### END USER INPUT
import os
import numpy as np
import matplotlib.pyplot as plt
# create directory structure if not in place
if not os.path.isdir('../output/PLOTS'):
    os.mkdir('../output/PLOTS')
# create directory structure if not in place
if not os.path.isdir('../output/PLOTS/GRESH_COMPARE'):
    os.mkdir('../output/PLOTS/GRESH_COMPARE')
# ring name and location data
lab = ['6','5','4',r'$\alpha$',r'$\beta$',r'$\eta$',r'$\gamma$',
        r'$\delta$',r'$\varepsilon$']
rings_km = {'E':[41794.93, 42155.17, 42555.67, 44686.59,
       45673.33 , 47176.51, 47628.26, 48297.35, 51342.23],
       'I':[41871.35, 42304.04, 42556.28 , 44736.75,
       45640.74, 47176.28, 47621.59, 48300.57, 50796.85]}
width = {'E':[1.72,2.62,2.67,4.22,11.19,1.53,1.63,2.7,74.93],
        'I':[1.52,2.75,1.95,10.59,7.03,1.54,3.83,6.7,22.43]}
# Gresh 1989 ring IDs
rnames = ['6','5','4','A','B','N','G','D','E']

for dir in ['I','E']:
    for i,name in enumerate(rnames):

        # load in GRESH data
        gfile = '../data/URINGS_GRESH1989_050M/RU1P2X'+name+dir+'.TAB'
        rho_dg,tau_dg = np.loadtxt(gfile,usecols=(0,1),delimiter=',').T
        # load in
        refdir = '../output/'+dir+'/'+name+'/'
        rfile = 'VGR2_X43_'+dir+'_URING_'+name+'_TAU_'+res+'M_'+date+'_'+sn+'.TAB'
        rho,tau,thresh = np.loadtxt(refdir+rfile,usecols=(0,6,8),delimiter=',').T

        plt.plot(rho_dg,tau_dg,'-r',label='Gresh 1989 50 m',lw=1)
        plt.plot(rho,tau,'-k',lw=1,label='TC2019 '+res[-2:]+'m')
        plt.plot(rho,thresh,'--c',label=r'$\tau_{thres}$ at 50 m')

        r = [(rho>rings_km[dir][i]-width[dir][i])&(rho<rings_km[dir][i]+width[dir][i])]
        maxt = [np.max(tau[r])]
        r = [(rho_dg>rings_km[dir][i]-width[dir][i])&(rho_dg<rings_km[dir][i]+width[dir][i])]
        maxt += [np.max(tau_dg[r])]
        plt.ylim(np.max(maxt)+0.1,-0.5)
        if 2*width[dir][i] < 10 :
            plt.xlim(rings_km[dir][i]-5,rings_km[dir][i]+5)
        elif 2*width[dir][i] > rho_dg[-1]-rho_dg[0] :
            plt.xlim(rho_dg[0],rho_dg[-1])
        else:
            plt.xlim(rings_km[dir][i]-width[dir][i],rings_km[dir][i]+width[dir][i])
        plt.xlabel(r'$\rho$')
        plt.ylabel(r'$\tau$')
        plt.title('Comparison of Gresh 1989 and TC Processing\nfor Uranus ring '+name+' during '+dir)
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig('../output/PLOTS/GRESH_COMPARE/VGR2_'+dir+'_URING_'+name+'_COMPARE_GRESH1989_'+res+'M.PDF')
        plt.close()
