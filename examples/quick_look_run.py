"""

quick_look_run.py

Purpose: Example rss_ringoccs/ 'quick-look' process script.

Revisions:
    2018 Sep 19 - jfong - original
"""

import sys
sys.path.append('../')
import rss_ringoccs as rss
sys.path.remove('../')
import pdb
import matplotlib.pyplot as plt

# ***** Begin user input *****
data_dir = '../output/Rev007/E/Rev007E_RSS_2005_123_X43_E/'
geo_file = data_dir + ''
cal_file = data_dir + ''
dlp_file = data_dir + ''

verbose = True
res_km = 1.0
inversion_range = 'all'
res_list = [1.5, 1.2, 0.75, 0.50]
Maxwell_xrange = [87400., 87600.]

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
    tau_inst = rss.diffcorr.DiffractionCorrection(dlp_inst, res,
            rng=Maxwell_xrange, verbose=verbose)
    rho = tau_inst.rho_km_vals - 87515.
    tau = tau_inst.tau_vals
    ax = axes[n]
    ax.plot(rho, tau, label=label)
    ax.set_ylabel('$\\tau$')
    ax.grid(True)
    ax.legend()

axes[0].set_title('Rev007E X43 Maxwell Ringlet Optical Depth \nReconstruction Resolution Comparison', fontsize=13)
axes[-1].set_xlabel('$\\rho$ - 87515 (km)')
axes[-1].set_xlim([-45.,45.])

pdb.set_trace()

