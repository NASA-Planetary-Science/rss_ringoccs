import numpy as np
import sys
sys.path.append('../')
import rss_uringoccs.calibration.calc_tau_thresh as ctt
from rss_uringoccs.calibration.resample_IQ import resample_IQ
from rss_uringoccs.tools.write_output_files import construct_filepath
sys.path.remove('../')
import phase_unwrap_params as unwrap
import pipeline_params as pipe
import matplotlib.pyplot as plt
import os
import pdb
### suppress any annoying warning messages
import warnings as wrn # for error suppression
wrn.filterwarnings('ignore')
np.seterr(all='ignore')

# function to perform third-order phase steering
def correct_phase(spm,phase,dir,ring_name):
    # signal sample rates
    dt = (spm[-1]-spm[0])/float(len(spm))
    n = int(1/dt)
    ### Initial fit
    # unwrapping phase
    phase += np.pi
    phase = np.unwrap(phase)
    ### Set unwrapping correction parameters
    tst  = unwrap.params[dir][ring_name]['T_OET_START']
    ten  = unwrap.params[dir][ring_name]['T_OET_END']
    prng = unwrap.params[dir][ring_name]['PHASE_RANGE']
    incr = unwrap.params[dir][ring_name]['2PI_INCR']
    ### Correct the unwrapping
    for i in range(len(tst)):
        if incr[i] > 0:
            phase[(spm>=tst[i])&(spm<ten[i])&(phase<prng[i])] += incr[i]*2.*np.pi
        else:
            phase[(spm>=tst[i])&(spm<ten[i])&(phase>prng[i])] += incr[i]*2.*np.pi
    ### "fit" with smoothing boxcar convolution
    kernel = np.ones(n)*dt
    k = len(kernel)
    fit = np.convolve(phase,kernel,'same')
    # correct edges with linear extrapolation
    fit[:k] = np.polyval(np.polyfit(spm[k:int(3*k/2)],fit[k:int(3*k/2)],2),spm[:k])
    fit[-k:] = np.polyval(np.polyfit(spm[int(-3*k/2):-k],fit[int(-3*k/2):-k],2),spm[-k:])
    ### Exclude ring regions to preserve diffraction pattern
    for rn in pipe.rnames:
        # get ring center
        ring_cent = pipe.rings_spm[dir][rn]
        # reject rings outside data slice
        if ring_cent < spm[0] or ring_cent > spm[-1] :
            continue
        # find edges of possible diffraction pattern using ring width
        # to infer the breadth of the diffraction pattern assuming
        # the pattern is 1 sec distant from the ring edges
        twidth = pipe.ring_widths[dir][rn]/16  # ring width in seconds
        rmin = np.argmin(abs(spm-(ring_cent-twidth-2)))
        rmax = np.argmin(abs(spm-(ring_cent+twidth+2)))
        # set fitting range to one second before and after
        fitargs = [i for i in range(rmin,rmin+5000)]
        fitargs += [i for i in range(rmax-5000,rmax)]
        # compute polynomial fit to phase adjacent to ring region
        par = np.polyfit(spm[fitargs],fit[fitargs],2)
        for i in range(rmin,rmax):
            fit[i] = np.polyval(par,spm[i])
    ### cusp between eta and gamma in ingress
    tcusp = 81335.9952
    if tcusp in spm :
        # mask to select 1 sec of data on each side of the cusp
        m1 = [(spm>81335)&(spm<tcusp)]
        m2 = [(spm>tcusp)&(spm<81337)]
        # polynomial fits to each data slice
        p1 = np.polyfit(spm[m1],phase[m1],2)
        p2 = np.polyfit(spm[m2],phase[m2],2)
        # replace "fit" with approximations at cusp
        for i in np.argwhere((spm>81335)&(spm<81337)) :
            if spm[i] < tcusp :
                fit[i] = np.polyval(p1,spm[i])
            else:
                fit[i] = np.polyval(p2,spm[i])
    # subtract fit to "flatten"
    diff = phase-fit
    # re-wrap, just to be safe, so that no values exceed [-pi,pi]
    dwrap= (diff+np.pi)%(2.*np.pi)-np.pi
    # smooth for plotting
    dsmooth = np.convolve(dwrap,np.ones(100)/100.,'same')
    dsmooth[:int(n/2)] = np.nan
    dsmooth[-int(n/2):] = np.nan

    # output file names for plotting
    file1,dir1 = construct_filepath(rev_info,'PHASE_FIT')
    file2,dir2 = construct_filepath(rev_info,'PHASE_FIT_ZOOM')
    file3,dir3 = construct_filepath(rev_info,'PHASE_STEER')

    # plot for posterity
    plt.plot(spm,phase,color='0.5',lw=1)
    plt.plot(spm,fit,'-k',lw=1)
    plt.xlabel('SPM (sec)')
    plt.ylabel('Phase (rad)')
    plt.title(dir+' '+ring_name+' Unwrapped Phase Fit')
    plt.savefig(dir1[0]+file1[0]+'.PDF')
    plt.close()

    # time limits
    ring_cent = pipe.rings_spm[dir][ring_name]
    twidth = pipe.ring_widths[dir][ring_name]/16  # ring width in seconds
    rmin = np.argmin(abs(spm-(ring_cent-twidth-2)))
    rmax = np.argmin(abs(spm-(ring_cent+twidth+2)))
    tmin = round(spm[rmin-10000])
    tmax = round(spm[rmax+10000])
    clp = [(spm>=tmin)&(spm<=tmax)]

    # plot for posterity
    plt.plot(spm,dwrap,color='0.5',lw=1)
    plt.plot(spm,dsmooth,'-k',lw=1)
    plt.xlabel('SPM (sec)')
    plt.ylabel('Phase (rad)')
    plt.xlim(tmin,tmax)
    plt.title(dir+' '+ring_name+' Residual Re-wrapped Steered Phase')
    #plt.show()
    plt.savefig(dir3[0]+file3[0]+'.PDF')
    plt.close()

    # return steered phase
    return dwrap

# function to redo power normalization
def normalize(rho,power):
    # mask array to exclude rings to preserve diffraction pattern
    ring_excl = np.array([True for i in range(len(rho))],dtype='bool')
    for rn in pipe.rnames:
        # get ring center
        ring_cent = pipe.rings_km[dir][rn]
        # reject rings outside data slice
        if ring_cent < rho[0] or ring_cent > rho[-1] :
            continue
        # find edges of possible diffraction pattern using ring width
        # to infer the breadth of the diffraction pattern assuming
        # the pattern is 16 km distant from the ring edges
        twidth = pipe.ring_widths[dir][rn]  # ring width
        rmin = np.argmin(abs(rho-(ring_cent-twidth-16)))
        rmax = np.argmin(abs(rho-(ring_cent+twidth+16)))
        for j in range(rmin,rmax+1):
            ring_excl[j] = False
    # trim off the ends
    ring_excl[:50] = False
    ring_excl[-50:] = False
    ring_excl[(power>3e3)] = False
    # polynomial fit
    p = np.polyfit(rho[ring_excl],power[ring_excl],9)
    pnorm_fit = np.polyval(p,rho)

    # output file names for plotting
    file1,dir1 = construct_filepath(rev_info,'POWER_NORM')

    # plot for posterity
    plt.plot(rho,power,color='0.75',lw=1)
    plt.plot(rho[ring_excl],power[ring_excl],color='0.5',lw=1)
    plt.plot(rho,pnorm_fit,'-k',lw=1)
    plt.ylim(-100,np.max(power[ring_excl])+500)
    plt.xlabel(r'$\rho$ (km)')
    plt.ylabel('Power (arb.)')
    plt.title(dir+' '+ring_name+' Freespace Power Fit')
    plt.savefig(dir1[0]+file1[0]+'.PDF')
    plt.close()

    return pnorm_fit,p

# create diffraction-limited profiles
# for each profile direction
for dir in ['I','E']:
    # set day-of-year based on profile direction
    if dir == 'I':
        doy = 24
    elif dir == 'E':
        doy = 25
    # create diffraction-limited profiles for each ring
    for ring_name in pipe.rnames:

        # naming convention for profile direction and ring
        fstart = '../output/'+dir+'/'+ring_name+'/VGR2_X43_'+dir+'_URING_'+ring_name
        # geo file
        geofile = fstart+'_GEO_'+pipe.geo_date+'_'+pipe.geo_sn+'.TAB'
        # cal file
        calfile = fstart+'_CAL_'+pipe.cal_date+'_'+pipe.cal_sn+'.TAB'

        # check ring name and files, and proceed if files and ring name exist
        if os.path.exists(geofile) and os.path.exists(calfile):
            # Rev info as gleaned from Gresh et al. 1989 reg: the Voyager 2
            # occultation event by Uranus on Jan 24, 1986 observed by the 70 m
            # receiver at the DSN station 43 in Canberra, Australia
            rev_info = { "rsr_file": 'None', "band": 'X', "year": '1986',
                "doy": str(doy), "dsn": '43', "occ_dir": dir, "planetary_occ_flag": 'None',
                "rev_num": 'None', "prof_dir": dir, "ring": ring_name, "spacecraft": 'VGR2',
                "phase_unwrap_params":unwrap.params[dir][ring_name] }

            if pipe.verbose:
                print('\nCreating DLP for Voyager 2 '+dir+'\noccultation of '+ring_name+' ring')
            # ring center
            spm_cent = pipe.rings_spm[dir][ring_name]
            # import geometry results
            if pipe.verbose:
                print('Importing GEO file')
                print('\t',geofile)
            geo = np.loadtxt(geofile,delimiter=',').T
            rsort = np.argsort(geo[3])
            # import calibration results
            if pipe.verbose:
                print('Importing CAL file')
                print('\t',calfile)
            spm,f_offset,f_sky,f_offit,I,Q = np.loadtxt(calfile,delimiter=',').T

            # corrected complex signal
            IQ_c = I+1j*Q
            # signal phase
            phase = -1. * np.arctan2(np.imag(IQ_c),np.real(IQ_c))
            # 3rd order phase steering
            phase_corr = correct_phase(spm,phase,dir,ring_name)

            # interpolate geometry to calibration data to get tau_thresh
            rho_geo = np.interp(spm,geo[0],geo[3])
            rho_dot = np.interp(spm,geo[0],geo[8])
            B_rad = np.interp(spm,geo[0],geo[6])*np.pi/180.

            # resample I and Q to radius
            rho_resamp,IQ_resamp = resample_IQ(rho_geo,IQ_c,pipe.dr_km_desired)
            oet = np.interp(rho_resamp,geo[3][rsort],geo[0][rsort])
            rho = rho_resamp
            IQ = IQ_resamp

            # interp to uniform radial sampling
            rho_dot = np.interp(rho,geo[3][rsort],geo[8][rsort])
            B = np.interp(rho,geo[3][rsort],geo[6][rsort])*np.pi/180.
            oet = np.interp(rho,geo[3][rsort],geo[0][rsort])
            ret =  np.interp(rho,geo[3][rsort],geo[1][rsort])
            set =  np.interp(rho,geo[3][rsort],geo[2][rsort])
            phi_rl = np.interp(rho,geo[3][rsort],geo[4][rsort])*np.pi/180.
            phi_ora = np.interp(rho,geo[3][rsort],geo[5][rsort])*np.pi/180.
            F = np.interp(rho,geo[3][rsort],geo[10][rsort])
            D = np.interp(rho,geo[3][rsort],geo[7][rsort])
            # sky frequency and normalized power
            f_sky = np.interp(oet,spm,f_sky)
            # redo power normalization
            pnorm,params = normalize(rho,abs(IQ)**2)
            power = (abs(IQ)**2)/abs(pnorm)
            # normal optical depth, phase, and threshold opticla depth
            tau = -np.sin(B)*np.log(power)
            phase = np.interp(oet,spm,phase_corr)
            rsort2 = np.argsort(rho_geo)

            # get the threshold optical depth from spectrogram and geometry
            tt_inst = ctt.calc_tau_thresh(spm,IQ_c,
                            np.interp(spm,geo[0],geo[3]),
                            np.interp(spm,geo[0],geo[8]),
                            np.interp(spm,geo[0],geo[6])*np.pi/180.,
                            params,pipe.rings_spm[dir][ring_name],
                            pipe.ring_widths[dir][ring_name],
                            res_km=pipe.dr_km_desired)
            # resample
            tau_thresh = np.interp(rho,rho_geo[rsort2],tt_inst.tau_thresh[rsort2])

            # write out a custom DLP file (DOES NOT MATCH PDS FORMAT!!)
            # create output file name
            if pipe.dr_km_desired < 0.01 :
                dr_str = '00' + str(int(pipe.dr_km_desired*1e3))
            elif pipe.dr_km_desired < 0.1 :
                dr_str = '0' + str(int(pipe.dr_km_desired*1e3))
            else:
                dr_str = str(int(pipe.dr_km_desired*1e3))
            # get radial range of ring profile needed for 20 m reconstruction
            # pad for highest resolution reconstruction assuming W = 2 F^2/dR
            # where dR is 2*dr_km_desired and F is the Fresnel scale
            pad = np.max(geo[10])**2 / pipe.dr_km_desired
            # compute min and max radii needed for reconstruction
            if pipe.ring_widths[dir][ring_name] > 25 :
                rho_min = round(pipe.rings_km[dir][ring_name]-pipe.ring_widths[dir][ring_name]-pad-1)
                rho_max = round(pipe.rings_km[dir][ring_name]+pipe.ring_widths[dir][ring_name]+pad+1)
            else:
                rho_min = round(pipe.rings_km[dir][ring_name]-26-pad)
                rho_max = round(pipe.rings_km[dir][ring_name]+26+pad)
            # write diffraction-limited profile to file over radius range
            # necessary for reconstruction
            outfile,outdir = construct_filepath(rev_info,'DLP_'+dr_str+'M')
            if pipe.verbose:
                print('Exporting DLP file')
                print('\t'+outdir[0]+outfile[0]+'.TAB')
            # write out DLP, CAL, and GEO info to TAB file
            out = open(outdir[0]+outfile[0]+'.TAB','w')
            for j in range(len(tau_thresh)):
                # if within reconstruction range
                if rho[j] >= rho_min and rho[j] <= rho_max :
                    vals = [rho[j],0.0,0.0,phi_rl[j],phi_ora[j],power[j],tau[j],
                            phase[j],tau_thresh[j],oet[j],ret[j],set[j],B[j],
                            rho_dot[j],F[j],D[j],f_sky[j]]
                    out.write(','.join(str(round(val,8)).rjust(18) for val in vals)+'\n')
            out.close()
        # in file name not recognized, say so
        else:
            print('One of the required GEO and CAL files not found.'+
                '\n'+'Check \'pipeline_params.py\' input dates and serial'+
                '\n'+'numbers for GEO and CAL against output files'+
                '\n'+'in \'./output/'+dir+'/'+ring_name+'/'+'\' to ensure'+
                '\n'+'these files exist:'+'\n\t'+geofile+'\n\t'+calfile)
