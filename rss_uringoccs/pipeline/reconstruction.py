import os
import sys
sys.path.append('../')
from rss_uringoccs.tools.history import write_history_dict
from rss_uringoccs.tools.write_output_files import write_output_files
sys.path.remove('../')
sys.path.append('../../')
import rss_ringoccs.diffrec.diffraction_correction as DiffCorr
sys.path.remove('../../')
import pipeline_params as pipe
import numpy as np
import matplotlib.pyplot as plt

###
### Class to read in values from DLP and create a faux DLP object
#       Note: this does not include any rev info or other LBL
#       information for the rss_ringoccs write_output_files. Any
#       output will need to be done manually.
###
class dlp(object):
    def __init__(self,file):

        data = np.loadtxt(file,delimiter=',').T

        self.rho_km_vals = data[0]              #rho_km_desired
        self.t_oet_spm_vals = data[9]           #spm_desired
        self.p_norm_vals = data[5]              #p_norm_vals
        self.phase_rad_vals = data[7]           #phase_rad_vals

        self.B_rad_vals = data[12]              #B_rad_vals_interp
        self.D_km_vals = data[15]               #D_km_vals_interp
        self.F_km_vals = data[14]               #F_km_vals_interp
        self.f_sky_hz_vals = data[16]           #f_sky_hz_vals_interp
        self.phi_rad_vals = data[4]             #phi_ora_rad_vals_interp
        self.t_ret_spm_vals = data[10]          #t_ret_spm_vals_interp
        self.t_set_spm_vals = data[11]          #t_set_spm_vals_interp
        self.phi_rl_rad_vals = data[3]          #phi_rl_rad_vals_interp
        self.rho_dot_kms_vals = data[13]        #rho_dot_kms_vals_interp
        self.rho_corr_pole_km_vals = data[1]    #rho_corr_pole_km_vals
        self.rho_corr_timing_km_vals = data[2]  #rho_corr_timing_km_vals
        self.tau_vals = data[6]                 #normal optical depth
        self.raw_tau_threshold_vals = data[8]   #threshold optical depth

        vars = {"DLP Data": file}
        kwds = {}
        self.history = write_history_dict(vars,kwds,__file__)

        if file.split('/')[2] == 'I':
            doy = 24
            pdir = 'I'
        elif file.split('/')[2] == 'E':
            doy = 25
            pdir = 'E'
        self.rev_info = {
            "rsr_file": "None",
            "band": 'X',
            "year": "1986",
            "doy": str(doy),
            "dsn": '43',
            "occ_dir": "BOTH",
            "planetary_occ_flag": "None",
            "rev_num": "None",
            "prof_dir": file.split('/')[2],
            "ring": file.split('/')[3]
        }

# for each profile direction
for dir in ['I','E']:
    # for each ring
    for ring_name in pipe.rnames:

        # reconstruct DLP filepath
        dlp_dir = '../output/' + dir + '/' + ring_name + '/'
        dlp_file = ('VGR2_X43_' + dir + '_URING_' + ring_name +
                    '_DLP_' + pipe.dlp_res + 'M_' + pipe.dlp_date +
                    '_' + pipe.dlp_sn + '.TAB')

        # rebuild DLP instance
        if pipe.verbose:
            print('Importing DLP file')
            print('\t',dlp_dir+dlp_file)
        dlp_inst = dlp(dlp_dir+dlp_file)

        # get radial range of ring profile
        if pipe.ring_widths[dir][ring_name] > 25 :
            rho_min = pipe.rings_km[dir][ring_name]-pipe.ring_widths[dir][ring_name]
            rho_max = pipe.rings_km[dir][ring_name]+pipe.ring_widths[dir][ring_name]
        else:
            rho_min = pipe.rings_km[dir][ring_name]-25
            rho_max = pipe.rings_km[dir][ring_name]+25

        # do reconstruction
        tau_inst = DiffCorr.DiffractionCorrection(dlp_inst,pipe.res,rng=[rho_min,rho_max],
                    res_factor=1, psitype=pipe.psitype, wtype=pipe.window,
                    fwd=False, norm=False, bfac=True, write_file=False,
                    verbose=pipe.verbose, sigma=2e-12)

        # write out to file
        if pipe.write_file :
            write_output_files(tau_inst)
        print('\n')
