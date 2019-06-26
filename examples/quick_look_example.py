"""
quick_look_run.py

Purpose: Example rss_ringoccs/ 'quick-look' process script.

"""

import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import matplotlib.pyplot as plt
import time

# ***** Begin user input *****
data_dir = '../output/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/'
geo_file = data_dir + 'RSS_2005_123_X43_E_GEO_YYYYMMDD_0001.TAB'
cal_file = data_dir + 'RSS_2005_123_X43_E_CAL_YYYYMMDD_0001.TAB'
dlp_file = data_dir + 'RSS_2005_123_X43_E_DLP_0100M_YYYYMMDD_0001.TAB'

verbose = True
res_km = 1.0
inversion_range = [87400., 87600.]  # Maxwell range
res_list = [2.0, 1.5, 1.0, 0.50, 0.25]

plot_center_km = 87515. # radius to center plot at

# ***** End user input *****

dlp_inst = rss.tools.ExtractCSVData(geo_file, cal_file, dlp_file,
        verbose=verbose)

nres = len(res_list)
fig, axes = plt.subplots(nres, 1, figsize=(8.5,11), sharex=True)
plt.subplots_adjust(hspace=0)
for n in range(nres):
    res = res_list[n]
    label = str(res*1000.) + 'm'
    print('Reconstructing at ', label, '...')
    start = time.time()
    tau_inst = rss.diffrec.DiffractionCorrection(dlp_inst, res,
            rng=[6.5e4,1.5e5], res_factor=0.75, psitype='Fresnel4',
            wtype='kbmd20', fwd=False, norm=True, bfac=True,
            write_file=False, verbose=True)
    end = time.time()
    print(res,' ',end-start)
    rho = tau_inst.rho_km_vals - plot_center_km
    tau = tau_inst.tau_vals
    ax = axes[n]
    ax.plot(rho, tau, label=label)
    ax.set_ylabel('$\\tau$')
    ax.grid(True)
    ax.legend()
    ax.set_ylim([-0.5,4.0])

axes[0].set_title('Rev007E X43 Maxwell Ringlet Optical Depth \nReconstruction Resolution Comparison', fontsize=13)
axes[-1].set_xlabel('$\\rho$ - ' + str(plot_center_km) + ' (km)')
axes[-1].set_xlim([-45.,45.])
plt.show()

"""
Revisions:
"""
