import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep,splev
from scipy.signal import decimate
from scipy.integrate import simps
import spiceypy as spice
import sys
sys.path.append('../')
import rss_uringoccs.reader.vgr_uranus_reader as vur
from rss_uringoccs.occgeo import occgeo_uranus as occgeo
from rss_uringoccs.tools.write_output_files import construct_filepath
sys.path.remove('../')
import pipeline_params as pipe
### suppress any annoying warning messages
import warnings as wrn # for error suppression
wrn.filterwarnings('ignore')
np.seterr(all='ignore')
# phase-detrend IQ
def correct_IQ(spm_vals,IQ_m,f_spm,f_offset_fit):
    # Interpolate frequeny offset fit to 0.1 second spacing, since
    # this makes the integration later more accurate
    dt = 0.01
    npts = round((f_spm[-1] - f_spm[0]) / dt)
    f_spm_ip = f_spm[0] + dt * np.arange(npts)
    f_off_ip = np.interp(f_spm_ip,f_spm,f_offset_fit)

    # Integration of frequency offset fit to get phase detrending function.
    # Then interpolated to same SPM as I and Q
    f_detrend_ip = np.zeros(len(f_off_ip))
    for i in range(len(f_off_ip)):
        f_detrend_ip[i] = simps(f_off_ip[:i+1],x=f_spm_ip[:i+1])
    f_detrend_ip_rad = f_detrend_ip * (2.0 * np.pi)-np.pi
    f_detrend_rad = np.interp(spm_vals,f_spm_ip,f_detrend_ip_rad)

    # Apply detrending function
    IQ_c = IQ_m * np.exp(-1j * f_detrend_rad)

    return IQ_c

# ring reference frames
frames = ['6','5','4','ALPHA','BETA','ETA','GAMMA','DELTA','LAMBDA','EPSILON']

tot_st = time.time()
# for each direction
for dir in ['I','E']:
    # information relevant to profile direction
    if dir == 'I' :
        #profdir = 'ingress'
        fbreak = 81336.0
        doy = 24
    elif dir == 'E' :
        #profdir = 'egress'
        fbreak = 5709.0
        doy = 25

    # Rev info as gleaned from Gresh et al. 1989 reg: the Voyager 2
    # occultation event by Uranus on Jan 24, 1986 observed by the 70 m
    # receiver at the DSN station in Canberra, Australia
    rev_info = { "rsr_file": 'None', "band": 'X', "year": '1986',
        "doy": str(doy), "dsn": '43', "occ_dir": dir, "planetary_occ_flag": 'None',
        "rev_num": 'None', "prof_dir": dir, "ring": '' }

    # read in file
    stt = time.time()
    read_inst = vur.VGRUranusReader(pipe.rawpath,dir)
    spm_raw = read_inst.spm
    I = read_inst.I
    Q = read_inst.Q
    ent = time.time()
    if pipe.verbose:
        print('Read Time:'.rjust(24),round(ent-stt,3),' sec')

    # read in offset frequency and sky frequency
    f_spm,f_offset,f_recon_coarse = np.loadtxt('../tables/URING_'+dir+'_FOF.CSV',
                                delimiter=',').T
    # spline the results
    fof_coef = splrep(f_spm[(f_offset<-1e4)],f_offset[(f_offset<-1e4)])
    f_recon_coef = splrep(f_spm,f_recon_coarse)
    # fit frequency offset with 9th order polynomial
    stt = time.time()
    m1 = [(f_spm<fbreak)&(f_offset<-1e4)]
    p1 = np.polyfit(f_spm[m1],f_offset[m1],9)
    m2 = [(f_spm>fbreak)&(f_offset<-1e4)]
    p2 = np.polyfit(f_spm[m2],f_offset[m2],9)
    f_offset_fit = np.polyval(p1,spm_raw[(spm_raw<fbreak)])
    f_offset_coarse = np.polyval(p1,f_spm[(f_spm<fbreak)])
    f_offset_fit = np.append(f_offset_fit,np.polyval(p2,spm_raw[(spm_raw>=fbreak)]))
    f_offset_coarse = np.append(f_offset_coarse,np.polyval(p2,f_spm[(f_spm>=fbreak)]))
    # coarsely sampled sky frequency
    f_sky_coarse = f_recon_coarse + f_offset_coarse
    ent = time.time()
    if pipe.verbose:
        print('Freq Offset Fit Time:'.rjust(24),round(ent-stt,3),' sec')

    # process each ring
    for rind,ring_name in enumerate(pipe.rnames):

        # update rev info
        rev_info['ring'] = ring_name

        #### GEOMETRY
        ra_deg = 76.5969
        dec_deg = 15.1117
        nhat_p = spice.radrec(1., ra_deg*spice.rpd(), dec_deg*spice.rpd())
        f0 = 8420430456.1
        f_sky_hz_vals = f_recon_coarse#np.zeros(len(f_spm)) + f0
        geo_inst = occgeo.Geometry(1986,doy, 'DSS-43', f_sky_hz_vals, f_spm,
                'URANUS', 'VOYAGER 2', pipe.kernels, rev_info, pt_per_sec=1.,
                ref='B1950', nhat_p=nhat_p, ring_frame='URING_'+frames[rind],
                verbose=pipe.verbose, write_file=True)

        #####
        if pipe.verbose:
            print('Calibrating Voyager 2 '+dir+'\noccultation of '+ring_name+' ring')
        ### set spm limits
        spm0 = pipe.rings_spm[dir][ring_name]-100
        spm1 = pipe.rings_spm[dir][ring_name]+100
        window = [(spm_raw>=spm0)&(spm_raw<=spm1)]
        # grab data relevant to ring
        spm = spm_raw[window]
        IQ = I[window]+1j*Q[window]
        ring_ind = np.argmin(abs(spm-pipe.rings_spm[dir][ring_name]))
        # frequency offset fitting
        #f_offset_fit_c = np.polyval(np.polyfit(spm,splev(spm,fof_coef),2),spm)

        # phase-detrend  raw signal
        IQ_c = correct_IQ(spm,IQ,spm,f_offset_fit[window])#IQ_steer[window]
        spm_c = spm

        # decimate now, for ease of computation
        stt = time.time()
        decif = [5,2] # decimation factors
        d = 1.0 # total factor of decimation
        for f in decif:
            IQ_c = decimate(IQ_c,f,ftype='fir',zero_phase=True)
            d *= f
        # correct funky 1 sec offset in alpha egress data
        if ring_name == 'A' and dir == 'E':
            spm_c = spm[0]+1+np.arange(len(IQ_c))*d/5e4
        # otherwise, times are generated without 1 sec offset
        else:
            spm_c = spm[0]+np.arange(len(IQ_c))*d/5e4
        # find ring center in time series
        ring_ind_c = np.argmin(abs(spm_c-pipe.rings_spm[dir][ring_name]))
        # spline to sample offset frequency at the same rate as IQ_c
        f_offset_c = splev(spm_c,fof_coef)
        f_offset_fit_c = np.interp(spm_c,spm,f_offset_fit[window])
        f_recon = splev(spm_c,f_recon_coef)
        #phase = -1. * np.arctan2(np.imag(IQ_c),np.real(IQ_c))
        ent = time.time()
        # output new sample rate
        if pipe.verbose:
            print('Resampling Time:'.rjust(24),round(ent-stt,3),' sec')
            print('Calibration Sampling:'.rjust(24),5e4/d,' Hz')
        # compute sky frequency
        f_sky = f_recon + f_offset_fit_c

        # write to custom CAL file
        stt = time.time()
        #  create output file
        outfile,outdir = construct_filepath(rev_info,'CAL')
        if pipe.verbose:
            print('Output file:',outdir[0]+outfile[0]+'.TAB')
        out = open(outdir[0]+outfile[0]+'.TAB','w')
        for i in range(len(spm_c)):
            I_c = np.real(IQ_c[i])
            Q_c = np.imag(IQ_c[i])
            row = [spm_c[i],f_offset_c[i],f_sky[i],f_offset_fit_c[i],I_c,Q_c]
            out.write(', '.join(str(round(val,8)).rjust(18) for val in row)+'\n')
        out.close()
        ent = time.time()
        if pipe.verbose:
            print('Writeout Time:'.rjust(24),round(ent-stt,3),' sec\n')

    print('\n')
tot_en = time.time()
print('Total processing time:'.rjust(24),round(tot_en-tot_st,3)/60,' min')
