import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep,splev
from scipy.signal import decimate
from scipy.integrate import simps
import spiceypy as spice
import sys
sys.path.append('../')
import rss_ringoccs.reader.vgr_uranus_reader as vur
from rss_ringoccs.occgeo import occgeo_uranus as occgeo
from rss_ringoccs.tools.write_output_files import construct_filepath
from rss_ringoccs import scatter
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
                verbose=pipe.verbose, write_file=False)

        #####
        if pipe.verbose:
            print('Calibrating Voyager 2 '+dir+'\noccultation of '+ring_name+' ring')
        ### set spm limits
        spm0 = pipe.rings_spm[dir][ring_name]-500
        spm1 = pipe.rings_spm[dir][ring_name]+500
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

        # compute power
        power = abs(IQ_c)**2.

        # mask array to exclude rings to preserve diffraction pattern
        ring_excl = np.array([True for i in range(len(spm_c))],dtype='bool')
        for rn in pipe.rnames:
            # get ring center
            ring_cent = pipe.rings_spm[dir][rn]
            # reject rings outside data slice
            if ring_cent < spm_c[0] or ring_cent > spm_c[-1] :
                continue
            # find edges of possible diffraction pattern using ring width
            # to infer the breadth of the diffraction pattern assuming
            # the pattern is 16 km distant from the ring edges
            twidth = pipe.ring_widths[dir][rn]/16  # ring width
            rmin = np.argmin(abs(spm_c-(ring_cent-twidth-1)))
            rmax = np.argmin(abs(spm_c-(ring_cent+twidth+1)))
            for j in range(rmin,rmax+1):
                ring_excl[j] = False
        # trim off the ends
        ring_excl[:50] = False
        ring_excl[-50:] = False
        ring_excl[(power>3e3)] = False
        # polynomial fit
        p = np.polyfit(spm_c[ring_excl],power[ring_excl],9)
        pnorm_fit = np.polyval(p,spm_c)

        # Get limits on data from the geometry
        spm_geo = geo_inst.t_oet_spm_vals
        rho_geo = geo_inst.rho_km_vals
        rho = np.interp(spm_c,spm_geo,rho_geo)
        rho_limits = [np.min(rho),np.max(rho)]
        scat_inst = scatter.Scatter(spm,spm_geo,rho_geo,spm_c,IQ_c,
                                    pnorm_fit,rev_info,
                                    rho_limits=rho_limits,
                                    nstack=int(32))

    print('\n')
