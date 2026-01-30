# rss_ringoccs_local_tools.py

# Revisions:
#   2026 Jan 28- rfrench 

# def append_file_to_file(source_file_path, destination_file_path):
# def b_MTR33(f0=8427222034.34050,band='X',allan_dev_1sec=2.e-13,rhodot=10,W=100.):
# def check_substring_in_list(substring, list_):
# def compute_tau_threshold(cal_file,tau_inst,verbose=False):
# def data_from_inst(dlp_inst,verbose=False):
# def DeltaR_phi_MTR32(b,DeltaR_W=1.):
# def DeltaR_inf_MTR3(b,Weff,F):
# def demo_e2e_event(Rev='Rev007',direction='E',DSN='X43',dr_km_desired=0.25,res_km=1.0,res_factor=0.75,
# def demo_e2e_event_loop(Rev='Rev007',direction='E',DSN='X43',dr_km_desired=0.25,
# def demo_Rev007_X43E_Maxwell(show=False,_16kHz=False,
# def demo_Rev007_X43E_Maxwell_loop(show=False,_16kHz=False,
# def demo_Rev007_X43I_Maxwell(show=False,_16kHz=False,
# def demo_Fring(Rev = 'Rev054',direction = 'I',DSN='K55',show=False,dr_km_desired=0.075,psitype='fresnel',
# def demo_StrangeRinglet(Rev = 'Rev067',direction = 'E',DSN='X14',
# def demo_Rev133I_X25_Cripples_loop(show=False,
# def demo_Rev125I_X34_Cripples_loop(show=False,
# def demo_Rev007E_X43_W7493_loop(show=False,_16kHz=False,
# def demo_Rev133I_X25_W7494_loop(show=False,_16kHz=False,
# def differential_opacity_plot(rev ='007',ID='Rev007E_rgf',direc ='E',
# def dtau_int(a,dQ,q=3.1):
#def event_name_from_rev_info(rev_info): 
#def event_name_from_taufilepath(taufilepath):    
# def featurelist_to_dict(infile = global_path_to_local_tables + 'my_feature_list_short_Cring.csv'):
# def find_first_nan_index(arr):
# def find_index(my_list, my_string,verbose=False):
# def flatten_extend(matrix):
# def format_time(name,rev,band,dsn,direc,dlp_res_used,psitype,res_km,tstart,tend,MINUTES=False):
# def get_all_RSR_files(local_path_to_tables=global_path_to_tables,force=False,download=True,silent=False):
# def get_all_SC_kernels(kernels_list = global_path_to_tables+'e2e_kernels.ker',
# def get_CORSS_8001_file(CORSS_8001_filepath,local_path_to_data=global_path_to_data, 
# def get_CORSS_8001_TABfiles(CORSS_8001_all_filepaths=None,Rev=None,\
# def get_CORSS_8001_TAUfile(rev_info,local_path_to_data = global_path_to_data):
# def get_CORSS_8001_XKa_TABfiles(CORSS_8001_all_filepaths=None,Rev=None,\
# def get_dBHz(_16kHz=True,Rev='007',direction='I',DSN='X43',
# def get_dBHz_from_rsr_file(rsr_file,
# def get_files_from_web(files,local_path,webURL,force=False,silent=True):
# def get_kernels_from_web(kernels,local_path_to_kernels=global_path_to_kernels,force=False,silent=False):
# get_loaded_kernels(basename=True)
# def get_processor_info():
# def get_psitype(labelfile):
# def get_rev_info_from_tau(taufile,rsr_file='"UNKNOWN"'):
# def get_RSRfiles_from_web(RSRfiles,local_path_to_data=global_path_to_data,force=False,silent=False):
# def get_trajectory_correction_coefficients(rev_info,fit_number=1,verbose=False,
# def is_res_km_valid(dlp_file,res_factor,res_km):
# def kernels_for_demo(verbose=False): 
# def kernel_is_loaded(kernel):
# def load_kernels_for_taufile(taufile,silent=True,load_absent_only=True): 
# def merge_sorted_tabfiles(indir,psitype=None,tmp2output=True,
# def min_required_dlp_res_km(res_km,res_factor):
# def min_required_dr_km_desired(res_km,resolution_factor):
# def min_valid_res_km(dlp_file,res_factor):
# def MTR86_Fig10(figsize=(6,7),figfile = global_path_to_local_figs+'MTR86_Fig10.jpg'):
# def MTR86_Fig11(figsize=(6,6),figfile = global_path_to_local_figs+'MTR86_Fig11.jpg'):
# def MTR86_Fig11_Cassini(figsize=(6,6),title='MTR86 Fig. 11 Cassini Rev137E',
# def n(a,q):
# def particle_size_models(w_K = 9e6,w_X = 3.6e7,w_S = 13.e7,q_vals = [2.8,3.0,3.2,3.4],
# def pick_SC_kernel(CORSS_TABfile,kernels_list = global_path_to_tables+'e2e_kernels.ker'):
# def pick_RSR_file(_16kHz=False,Rev='Rev007',direction='E',DSN='X43',
# def plot_comparisons(program,results,show=False,y_offset=0.1,dy_offset=0.,plot_phase=False,tau=None,
# def printnow(CRbefore=True,CRafter=True):
# def print_featurelist(d):
# def radius_correction_pole(tau_file,fit_number=1,verbose=False):
# def radius_correction_pole_from_tau_inst(dlp_file,tau_inst,fit_number=1,verbose=False):
# def radius_correction_trajectory(tau_file,fit_number=1,verbose=False):
# def radius_correction_trajectory_from_tau_inst(dlp_file,tau_inst,fit_number=1,verbose=False):
# def rkm_tau_from_TAUfile(path_to_tau_file)::
# def rkm_tau_tanBeff_from_TAUfile(path_to_tau_file):
# def rkm_tau_drs_tanBeff_from_TAUfile(path_to_tau_file):
# def RSRfiles_for_demo(_16kHz=True): 
# def runloop(figfile=None,PRINT_TIME=False,include_history=False,verbose=False,xlim_run=(None,None),
# def set_tau_history(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
# def set_tau_history_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
# def spm2date(year,doy,spmval):
# def subdivide_inversion_range(inversion_range, drange):
# def tau_from_CSV(geo_file,cal_file,dlp_file,Rev,direction,DSN,res_km=1.0,res_factor=0.75,   
# def tau_from_tau_inst(tau_inst):
# def tau_int(a,Q,q=3.1):
# def trim_dlp_file(dlp_file,processing_range,path_to_output=global_path_to_local_output,
# def update_file_version(pathtofile):
# def update_taufile_radius_corrections(taufile,update_dr_pole=True,update_dr_trajectory=True,
# def write_tau_files(tau_inst,geo_inst,cal_inst,dlp_inst,dlp_file,tstart,tend,
# def write_tau_files_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                              
import rss_ringoccs_global_paths # required to allow program to run from other directories

global global_path_to_rss_ringoccs
global_path_to_rss_ringoccs = rss_ringoccs_global_paths.global_path_to_rss_ringoccs

global global_path_to_local # note shortened name in next line
global_path_to_local = rss_ringoccs_global_paths.global_path_to_rss_ringoccs_local

# auxiliary path definitions based on above two globals

global global_path_to_data
global_path_to_data = global_path_to_rss_ringoccs + 'data/'

global global_path_to_output
global_path_to_output = global_path_to_rss_ringoccs + 'output/'

global global_path_to_tables
global_path_to_tables = global_path_to_rss_ringoccs + 'tables/'

global global_path_to_kernels
global_path_to_kernels = global_path_to_rss_ringoccs + 'kernels/'

global global_path_to_demo
global_path_to_demo = global_path_to_rss_ringoccs + 'demo/'

global global_path_to_demo_figs
global_path_to_demo_figs = global_path_to_rss_ringoccs + 'demo/figs/'

global global_path_to_local_output
global_path_to_local_output = global_path_to_local + 'output/'

global global_path_to_local_data
global_path_to_local_data = global_path_to_local + 'data/'

global global_path_to_local_figs
global_path_to_local_figs = global_path_to_local + 'program/figs/'

global global_path_to_local_picklefiles
global_path_to_local_picklefiles = global_path_to_local + 'picklefiles/'

global global_path_to_local_tables
global_path_to_local_tables = global_path_to_local + 'tables/'

global global_path_to_local_tmp
global_path_to_local_tmp = global_path_to_local + 'tmp/'

import rss_ringoccs
import rss_ringoccs as rss

from rss_ringoccs.tools import write_output_files
from rss_ringoccs.tools.pds3_tau_series import get_rev_info_from_dlp,set_tau_history_from_csv,add_psitype_to_taufiles

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
import datetime
import fcntl
import glob
import math
import matplotlib.patches as pch
import matplotlib.pyplot as plt
from multiprocessing import Process
from multiprocessing import cpu_count
import numpy as np
import os
import platform
import re
import scipy.signal
from scipy import interpolate
from scipy.integrate import simpson
import shutil
import spiceypy as spice
import subprocess
import sys
import time
import traceback
import types
from wcmatch import glob

def append_file_to_file(source_file_path, destination_file_path):
    """Appends the content of the source file to the destination file."""
    try:
        with open(source_file_path, 'r') as source_file:
            source_content = source_file.read()
        with open(destination_file_path, 'a') as destination_file:
            destination_file.write(source_content)
    except FileNotFoundError:
        print(f"Error: One or both files not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
def b_MTR33(f0=8427222034.34050,band='X',allan_dev_1sec=2.e-13,rhodot=10,W=100.):
    '''
    MTR86 eq 33
    '''
    omega0 = 2*np.pi*f0
    Weff = W/1.65 # by MTR86 p. 132
    b = (omega0 * allan_dev_1sec)**2 * Weff/(2*rhodot)
    return b
def calc_rho_km(et_vals, planet, spacecraft, dsn, kernels=None,
        ring_frame=None):
    """
    Calculate the distance between Saturn center to ring intercept point.

    Arguments
        :et_vals (*np.ndarray*): Array of observed event times in ET sec.
        :planet (*str*): Planet name
        :spacecraft (*str*): Spacecraft name
        :dsn (*str*): DSN observing station ID

    Keyword Arguments
        :kernels (*str* or *list*): List of NAIF kernels, including path.
        :ring_frame (*str*): Ring plane frame. Default is the equatorial
                             frame, (e.g. 'IAU_SATURN')

    Returns
        :rho_km_vals (*np.ndarray*): Array of ring intercept points in km.
    """

    # Calculate rho vector
    rho_vec_km, t_ret_et = calc_rho_vec_km(et_vals, planet, spacecraft, dsn,
            kernels=kernels, ring_frame=ring_frame)

    # Compute magnitude of each vector
    rho_km_vals = [spice.vnorm(vec) for vec in rho_vec_km]

    return np.asarray(rho_km_vals)

def calc_rho_vec_km(et_vals, planet, spacecraft, dsn, ref='J2000', kernels=None,
        verbose=False, ring_frame=None):
    """
    This calculates the position vector of the ring intercept point from the
    planet center in J2000 frame.

    Arguments
        :et_vals (*np.ndarray*): Array of earth-received times in ET sec
        :planet (*str*): Name of planet
        :spacecraft (*str*): Name of spacecraft
        :dsn (*str*): DSN observing station ID

    Keyword Arguments
        :kernels (*str* or *list*): Path to NAIF kernels
        :verbose (*bool*): Option for printing processing steps
        :ring_frame (*str*): Ring plane frame. Default is the equatorial
                             frame, (e.g. 'IAU_SATURN')
        :ref (*str*): Reference frame to be used in spiceypy calls. Default
                      is 'J2000'

    Returns
        :rho_vec_km_vals (*list*): List of 3xN np.ndarrays of the planet
            center to ring intercept point position vector in J2000 frame
        :t_ret_et_vals (*np.ndarray*): Array of ring event times in ET seconds.

    References
        #. Ring intercept point calculation using a dynamical frame.
           See [NAIF]_ page 19.
    """

    if kernels:
        spice.kclear()
        spice.furnsh(kernels)

    planet_id = spice.bodn2c(planet)

    npts = len(et_vals)
    rho_vec_km_vals = []
    t_ret_et_vals = np.zeros(npts)
    
    # Replace Saturn radii in kernel pool with values that represent a
    #   flat ellipsoid, which we will use as the ring plane
    body = planet_id
    item = 'RADII'
    dim = 3
    radii = spice.bodvcd(body, item, dim)

    new_radii = [1.e6, 1.e6, 1.e-5]

    # Construct naif radii code in form: 'BODY699_RADII'
    
    planet_naif_radii = 'BODY'+str(planet_id)+'_RADII'

    name = planet_naif_radii
    dvals = new_radii
    spice.pdpool(name, dvals)
    iau_planet = 'IAU_'+planet.upper()
    if ring_frame is None:
        frame = iau_planet

    else:
        frame = ring_frame

    if verbose:
        print('frame',frame)

    for n in range(npts):
        et = et_vals[n]

        # Compute spacecraft position relative to dsn
        targ = spacecraft
        #ref = 'J2000'
        abcorr = 'CN'
        obs = dsn
        starg, ltime, = spice.spkpos(targ, et, ref, abcorr, obs)

        nhat_sc2dsn = spice.vhat(starg)

        # Compute intersection of vector with ring plane and time epoch
        #   of intersection in ET secs (ring event time)

        method = 'Ellipsoid'
        target = planet
        # fixref = iau_planet
        fixref = frame
    #    fixref = iau_planet
        abcorr = 'CN'
        obsrvr = dsn
        #dref = 'J2000'
        dref = ref
        dvec = nhat_sc2dsn
        #print('about to call sincpt...fixref=',fixref)
        spoint, trgepc, srfvec = spice.sincpt(method, target, et, fixref, abcorr, obsrvr, dref, dvec)

        t_ret_et_vals[n] = trgepc

        # Convert ring plane intercept to J2000 frame
        frame_from = fixref
        frame_to = dref
        etfrom = trgepc
        etto = et
        xform = spice.pxfrm2(frame_from, frame_to, etfrom, etto)

        rho_vec_km = spice.mxv(xform, spoint)
        rho_vec_km_vals.append(rho_vec_km)

    # Restore old original values of RADII to kernel pool
    name = planet_naif_radii
    dvals = radii[1].tolist()
    spice.pdpool(name, dvals)

    return rho_vec_km_vals, t_ret_et_vals

def check_substring_in_list(substring, list_):
  """Checks if a substring is present in any element of a list.

  Args:
    substring: The substring to search for.
    list_: The list of strings to check.

  Returns:
    True if the substring is found in any element, False otherwise.
  """
  for element in list_:
    if substring in element:
      return True
  return False
    
def compute_tau_threshold(cal_file,tau_inst,verbose=False):
    '''
    Compute tau threshold using independently calculated dBHz and other information from PDS files.
    Plot the results if requested.
    '''
    CALfile = cal_file
    CALfileLBL = CALfile.replace('TAB','LBL')
    f = open(CALfileLBL, "r")
    contents = f.read().splitlines() # strips \n at end of each line
    f.close()
    rsr_file = 'NOT FOUND'
    for line in contents:
        if "SOURCE_PRODUCT_ID" in line:
            rsr_file = line[len('SOURCE_PRODUCT_ID               = "'):-1]
            if verbose:
                print(line)
                print("rsr_file:",rsr_file)
            break
    
    if verbose:
        print(cal_file)
        print(rsr_file)
        
    dCAL = np.genfromtxt(CALfile,delimiter=',')
    oetCAL = dCAL[:,0]
    pwrCAL = dCAL[:,3]
    pwrCAL /= np.nanmax(pwrCAL)

    dBHz = get_dBHz_from_rsr_file(rsr_file)
    if verbose:
        print('got dBHz = ',dBHz)
    
    pwrCAL *= 10.**(dBHz/10)
    SNR0 = pwrCAL    
    
    rTAU = tau_inst.rho_km_vals 
    tauTAU = tau_from_tau_inst(tau_inst) 
    oetTAU = tau_inst.t_oet_spm_vals 
    BdegTAU = tau_inst.B_deg_vals
    resolution = tau_inst.input_resolution_km
    resolution_factor = tau_inst.resolution_factor
    mu = np.abs(np.sin(np.radians(BdegTAU[0])))
    rhodotTAU = tau_inst.rho_dot_kms_vals

    if rhodotTAU[0] <0 :
        frCAL = interpolate.interp1d(oetTAU[::-1],rTAU[::-1],bounds_error=False,fill_value='extrapolate')
        rCAL = frCAL(oetCAL[::-1])
        rCAL = rCAL[::-1]
        # plt.plot(rTAU[::-1],oetTAU[::-1])
        # plt.plot(rCAL,oetCAL)
        # plt.show()
    else:
        frCAL = interpolate.interp1d(oetTAU,rTAU,bounds_error=False,fill_value='extrapolate')
        rCAL = frCAL(oetCAL)
    # print('rCAL:',rCAL)
    # print('oetCAL:',oetCAL)
    rhodotCAL = np.gradient(rCAL,oetCAL,edge_order=2)

    # plt.plot(rTAU,rhodotTAU,label='TAU')
    # plt.plot(rCAL,rhodotCAL,label='CAL')
    # plt.legend()
    # plt.show()
    DeltaRW = resolution * resolution_factor# processing resolution
    #tauth_calc = -mu*np.log(1.205*0.32*1.65/(pwrCAL) * np.abs(rhodotCAL) / DeltaRW)
    
    tauth_calc = mu * np.log(1.571733936879165 * SNR0 * DeltaRW / np.abs(rhodotCAL)) 
    # print('tauth_calc:',tauth_calc)
    # print('rhodotCAL:',rhodotCAL)

    ftauthCAL = interpolate.interp1d(rCAL,tauth_calc,bounds_error=False,fill_value='extrapolate')
    tauth_compare = ftauthCAL(rTAU)

    # plt.plot(rTAU,tauth_compare)
    # plt.title('threshold optical depth')
    # plt.show()

    return tauth_compare
    
def compute_PDS_tau_threshold(Rev='Rev007',direction='E',DSN='X43',tau_match='TAU_01KM',resolution=1,
            xlim=(None,None),ylim=(None,None),res_factor=0.75,figsize=(8,6),plot=True,verbose=False):
    '''
    Compute tau threshold using independently calculated dBHz and other information from PDS files.
    Plot the results if requested.
    '''
    # make sure we have TAU, CAL LBL and TAB files 
    CORSS_8001_TAUfile = get_CORSS_8001_TABfiles(Rev=Rev,DSN=DSN,direction=direction)
    if verbose:
        print('CORSS_8001_TAUfile',CORSS_8001_TAUfile)
    local_CORSS_8001_TAUfile = global_path_to_data + CORSS_8001_TAUfile
    CORSS_8001_TAU_LBL = get_CORSS_8001_file(CORSS_8001_TAUfile.replace('TAB','LBL'),silent=True)
    
    temp = CORSS_8001_TAUfile.replace(tau_match,'CAL')
    temp = temp.replace('TAB','LBL')
    CORSS_8001_CALfile = get_CORSS_8001_file(CORSS_8001_TAUfile.replace(tau_match,'CAL'),silent=True)
    CORSS_8001_CAL_LBL = get_CORSS_8001_file(temp,silent=True)
    
    dTAU = np.genfromtxt(local_CORSS_8001_TAUfile,delimiter=',')
    rTAU = dTAU[:,0] 
    tauTAU = dTAU[:,6] 
    tauthTAU = dTAU[:,8] 
    oetTAU = dTAU[:,9] 
    retTAU = dTAU[:,10] 
    BdegTAU = dTAU[:,12]
    mu = np.abs(np.sin(np.radians(BdegTAU[0])))
    rhodotTAU = np.gradient(rTAU,oetTAU,edge_order=2) # ignore difference from retTAU
    
    if verbose:
        print(CORSS_8001_CALfile)
    dCAL = np.genfromtxt(CORSS_8001_CALfile,delimiter=',')
    oetCAL = dCAL[:,0]
    pwrCAL = dCAL[:,3]
    pwrCAL /= np.nanmax(pwrCAL)
    corr34M = 20.*np.log10(70./34)
    dBHz = get_dBHz(Rev=Rev,direction=direction,DSN=DSN)
    pwrCAL *= 10.**(dBHz/10)
    SNR0 = pwrCAL
    if direction == 'I':
        frCAL = interpolate.interp1d(oetTAU[::-1],rTAU[::-1],bounds_error=False,fill_value='extrapolate')
        rCAL = frCAL(oetCAL[::-1])
        rCAL = rCAL[::-1]
        # plt.plot(rTAU[::-1],oetTAU[::-1])
        # plt.plot(rCAL,oetCAL)
        # plt.show()
    else:
        frCAL = interpolate.interp1d(oetTAU,rTAU,bounds_error=False,fill_value='extrapolate')
        rCAL = frCAL(oetCAL)
    # print('rCAL:',rCAL)
    # print('oetCAL:',oetCAL)
    rhodotCAL = np.gradient(rCAL,oetCAL,edge_order=2)

    # plt.plot(rTAU,rhodotTAU,label='TAU')
    # plt.plot(rCAL,rhodotCAL,label='CAL')
    # plt.legend()
    # plt.show()
    DeltaRW = resolution * res_factor# processing resolution
    #tauth_calc = -mu*np.log(1.205*0.32*1.65/(pwrCAL) * np.abs(rhodotCAL) / DeltaRW)
    
    tauth_calc = mu * np.log(1.571733936879165 * SNR0 * DeltaRW / np.abs(rhodotCAL)) 
    # print('tauth_calc:',tauth_calc)
    # print('rhodotCAL:',rhodotCAL)

    ftauthCAL = interpolate.interp1d(rCAL,tauth_calc,bounds_error=False,fill_value='extrapolate')
    tauth_compare = ftauthCAL(rTAU)
    if plot:
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(rTAU,tauthTAU,label='PDS '+Rev+direction+'_'+DSN+' '+str(resolution)+ ' km')
        ax.plot(rTAU,tauth_compare,label='comparison')
        if ylim == (None,None):
            ylim = (np.min(tauthTAU)-.5,np.max(tauthTAU)+.5)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        ax.legend()
        ax.set_title(os.path.basename(CORSS_8001_TAUfile))
        plt.show()
    return rTAU,tauthTAU,tauth_compare
    
def data_from_inst(dlp_inst,verbose=False):
    # if verbose:
    #     all_attributes = dir(geo_inst)
    # #    print(f"All attributes: {all_attributes}\n")
    #     geo_attributes = [attr for attr in all_attributes if not attr.startswith('__')]
    #     print("GEO attributes:")
    #     for at in geo_attributes:
    #         try:
    #             print(at,len(getattr(geo_inst,at)))
    #         except:
    #             print(at)
    
    #     all_attributes = dir(cal_inst)
    # #    print(f"All attributes: {all_attributes}\n")
    #     cal_attributes = [attr for attr in all_attributes if not attr.startswith('__')]
    #     print("CAL attributes:")
    #     for at in cal_attributes:
    #         try:
    #             print(at,len(getattr(cal_inst,at)))
    #         except:
    #             print(at)
    if verbose:
        all_attributes = dir(dlp_inst)
    #    print(f"All attributes: {all_attributes}\n")
        dlp_attributes = [attr for attr in all_attributes if not attr.startswith('__')]
        print("DLP attributes:")
        for at in dlp_attributes:
            try:
                print(at,len(getattr(dlp_inst,at)))
            except:
                print(at)

# we don't need GEO or CAL!

# GEO attributes:
# B_deg_vals 12787
# B_eff_deg_vals 12787
# D_km_vals 12787
# F_km_vals 12787
# R_imp_km_vals 12787
# _Geometry__get_naif_version
# _Geometry__remove_atmos_values
# _Geometry__remove_values_beyond_rings
# add_info
# atmos_occ_spm_vals 2877
# beta_vals 12787
# elev_deg_vals 12787
# freespace_km 16
# freespace_spm 16
# get_chord_ind
# get_profile_dir
# history 11
# ionos_occ_et_vals 2840
# ionos_occ_spm_vals 2840
# kernels 13
# naif_toolkit_version 9
# outfiles 1
# phi_ora_deg_vals 12787
# phi_rl_deg_vals 12787
# phi_rl_dot_kms_vals 12787
# rev_info 9
# rho_dot_kms_vals 12787
# rho_km_vals 12787
# rx_km_vals 12787
# ry_km_vals 12787
# rz_km_vals 12787
# split_ind
# t_oet_et_vals 12787
# t_oet_spm_vals 12787
# t_ret_spm_vals 12787
# t_set_et_vals 12787
# t_set_spm_vals 12787
# verbose
# verify_chord
# vx_kms_vals 12787
# vy_kms_vals 12787
# vz_kms_vals 12787

# CAL attributes:
# FORFIT_chi_squared
# FSPFIT_chi_squared
# IQ_c 12786000
# correct_IQ
# f_offset 841
# f_offset_fit_vals 12787
# f_sky_hz_vals 12787
# f_spm 841
# gaps 16
# history 11
# naif_toolkit_version 9
# outfiles 1
# p_free_vals 12787
# rev_info 9
# t_oet_spm_vals 12787

# DLP has
# DLP attributes:
# DLP attributes:
# B_rad_vals 20001
# D_km_vals 20001
# F_km_vals 20001
# _DiffractionLimitedProfile__interp_and_set_attr
# add_info
# create_dlps
# dr_km
# f_sky_hz_vals 20001
# history 11
# outfiles 1
# p_norm_vals 20001
# phase_rad_vals 20001
# phi_rad_vals 20001
# phi_rl_rad_vals 20001
# raw_tau_threshold_vals 20001
# rev_info 9
# rho_corr_pole_km_vals 20001
# rho_corr_timing_km_vals 20001
# rho_dot_kms_vals 20001
# rho_km_vals 20001
# rx_km_vals 20001
# ry_km_vals 20001
# rz_km_vals 20001
# t_oet_spm_vals 20001
# t_ret_spm_vals 20001
# t_set_spm_vals 20001
# tau_vals 20001

# data has attributes:
# B_deg_vals 20001
# D_km_vals 20001
# f_sky_hz_vals 20001
# history 9
# p_norm_vals 20001
# phase_deg_vals 20001
# phi_deg_vals 20001
# phi_rl_deg_vals 20001
# raw_tau_threshold_vals 20001
# raw_tau_vals 20001
# rev_info None
# rho_corr_pole_km_vals 20001
# rho_corr_timing_km_vals 20001
# rho_dot_kms_vals 20001
# rho_km_vals 20001
# rx_km_vals 20001
# ry_km_vals 20001
# rz_km_vals 20001
# t_oet_spm_vals 20001
# t_ret_spm_vals 20001
# t_set_spm_vals 20001
# tau_phase_deg_vals None
# tau_power_vals None
# tau_vals None

    class RSSdata:
        pass
    d = RSSdata() # populate this
    setattr(d,'B_deg_vals',dlp_inst.B_rad_vals*180/np.pi)

    setattr(d,'D_km_vals',dlp_inst.D_km_vals)
    setattr(d,'f_sky_hz_vals',dlp_inst.f_sky_hz_vals)
#     history = {'rss_ringoccs Version': '1.3', 'libtmpl Version': '0.1', 'Python Version': '3.14', 'Host Name': 'Loon', 'User Name': 'rfrench', 'Run Date': 'Thu Jan 15 00:16:18 2026', 'Operating System': 'macOS', 'Positional Args': {'geo': '/Volumes/dione_raid2/Research/PDART_RSS/rss_ringoccs/output/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO_20260115_0005.TAB', 'cal': '/Volumes/dione_raid2/Research/PDART_RSS/rss_ringoccs/output/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_CAL_20260115_0005.TAB', 'dlp': '/Volumes/dione_raid2/Research/PDART_RSS/rss_ringoccs/output/Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_DLP_0500M_20260115_0005.TAB'}, 'Keyword Args': {'tau': 'None', 'use_deprecated': False}}
# #    setattr(d,'history',history)
    setattr(d,'history',None)
    setattr(d,'p_norm_vals',dlp_inst.p_norm_vals)
    setattr(d,'phase_deg_vals',dlp_inst.phase_rad_vals*180/np.pi)
    setattr(d,'phi_deg_vals',dlp_inst.phi_rad_vals*180/np.pi)
    setattr(d,'phi_rl_deg_vals',dlp_inst.phi_rl_rad_vals*180/np.pi)
    setattr(d,'raw_tau_threshold_vals',dlp_inst.raw_tau_threshold_vals)
    raws_tau_vals = dlp_inst.raw_tau_threshold_vals*0
    setattr(d,'raw_tau_vals',raws_tau_vals)

    setattr(d,'rev_info',None)

    rho_corr_pole_km_vals = dlp_inst.raw_tau_threshold_vals*0
    setattr(d,'rho_corr_pole_km_vals',rho_corr_pole_km_vals)
    rho_corr_timing_km_vals = dlp_inst.raw_tau_threshold_vals*0
    setattr(d,'rho_corr_timing_km_vals',rho_corr_timing_km_vals)
    setattr(d,'rho_dot_kms_vals',dlp_inst.rho_dot_kms_vals)
    setattr(d,'rho_km_vals',dlp_inst.rho_km_vals)
    setattr(d,'rx_km_vals',dlp_inst.rx_km_vals)
    setattr(d,'ry_km_vals',dlp_inst.ry_km_vals)
    setattr(d,'rz_km_vals',dlp_inst.rz_km_vals)
    setattr(d,'t_oet_spm_vals',dlp_inst.t_oet_spm_vals)
    setattr(d,'t_ret_spm_vals',dlp_inst.t_ret_spm_vals)
    setattr(d,'t_set_spm_vals',dlp_inst.t_set_spm_vals)
    setattr(d,'tau_phase_deg_vals',None)
    setattr(d,'tau_power_vals',None)
    setattr(d,'tau_phase_deg_vals',None)
    setattr(d,'tau_vals',None)

    if verbose:
        all_attributes = dir(d)
        print(f"All attributes: {all_attributes}\n")
        for a in all_attributes:
            print(a)
        d_attributes = [attr for attr in all_attributes if not attr.startswith('__')]
        print("data attributes:")
        for at in d_attributes:
            try:
                print(at,len(getattr(d,at)))
            except:
                print(at)
    return d

def DeltaR_phi_MTR32(b,DeltaR_W=1.):
    '''
    MTR86 eq 32
    '''
    return DeltaR_W * (b**2/2)/(np.exp(-b) + b - 1)

def DeltaR_inf_MTR3(b,Weff,F):
    '''
    MTR86 eq 10 as a function of b
    '''
    return F**2*b/Weff

def demo_e2e_event(Rev='Rev007',direction='E',DSN='X43',dr_km_desired=0.25,res_km=1.0,res_factor=0.75,
    profile_range=[85000., 90000.],inversion_range=[87475,87560],psitype='fresnel',
    ring_frame='IAU_SATURN',ylim_tau=(None,None),
    wtype='kbmd20',include_CORSS_8001=True,plot=True,plot_phase=False,local_path_to_data =global_path_to_data,
    adjust_phase=True,y_offset=0,dy_offset=0,local_path_to_kernels = global_path_to_kernels, 
    local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output,
    local_path_to_demo_figs=global_path_to_demo_figs,save_figfile=True,_16kHz=False,
    decimate_16khz_to_1khz=False,decimate_16khz_to_2khz=False,
    decimate_16khz_to_4khz=False,decimate_16khz_to_8khz=False,
                decimate_50khz_to_25khz=False,
                decimate_50khz_to_10khz=False,
                decimate_50khz_to_5khz=False,
                decimate_50khz_to_1khz=False,
    verbose=True,silent=False,show=True,program='demo_e2e_event',use_CSVs=True,PerformInversion=True,
                  compute_tau_threshold_ = True, compute_dr_pole=True, compute_dr_trajectory=True,fit_number=1):
    '''
    perform end-to-end Diffraction reconstruction, beginning with RSR file, creating GEO,CAL,DLP, and TAU files
    '''
    if PerformInversion and (res_factor * res_km < 2*dr_km_desired):
            print("Unable to perform diffraction correction at requested resolution: res_factor * res_km < 2*dr_km_desired")
            print('res_factor',res_factor,'res_km',res_km,'dr_km_desired',dr_km_desired)
            print('Highest achievable resolution res_km = ',f'{2*dr_km_desired/res_factor:0.4f} km')
            return -1

    band = DSN[0]
    
    planet = 'Saturn' 
    spacecraft = 'Cassini'
    pnf_order = 3 # polynomial order for normalizing power vs time
    CORSS_8001_filepath =get_CORSS_8001_TABfiles(CORSS_8001_all_filepaths=None,Rev=Rev,\
                DSN=DSN,direction=direction,local_path_to_data=local_path_to_data,silent=silent)

    # set timer for job execution, included in output file header

    start_time = time.time()

    # get required kernel files if not already fetched

    kernels = get_kernels_from_web(kernels_for_demo(),\
                local_path_to_kernels=local_path_to_kernels,silent=silent)
    kernel_naif = pick_SC_kernel(CORSS_8001_filepath)
    kernelSC = get_kernels_from_web(kernel_naif,local_path_to_kernels=local_path_to_kernels,silent=True)[0]
    kernels.append(kernelSC)

    # get required data files for all available DSNs if not already fetched

    rsr_file = local_path_to_data + \
           pick_RSR_file(_16kHz=_16kHz,Rev=Rev,direction=direction,DSN=DSN,download=True,silent=True)

    # set up specifications for creating GEO/CAL/DLP files

    # Create instance with rsr file contents

    rsr_inst = rss.rsr_reader.RSRReader(rsr_file,
                                        decimate_16khz_to_8khz=decimate_16khz_to_8khz,
                                        decimate_16khz_to_4khz=decimate_16khz_to_4khz,
                                        decimate_16khz_to_2khz=decimate_16khz_to_2khz,
                                        decimate_16khz_to_1khz=decimate_16khz_to_1khz,
                                        decimate_50khz_to_25khz=decimate_50khz_to_25khz,
                                        decimate_50khz_to_10khz=decimate_50khz_to_10khz,
                                        decimate_50khz_to_5khz=decimate_50khz_to_5khz,
                                        decimate_50khz_to_1khz=decimate_50khz_to_1khz,
                                        verbose=True)

    # Create instance with geometry parameters

    geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft,
                kernels,
                ring_frame=ring_frame,
                local_path_to_tables=global_path_to_tables,
                local_path_to_output=global_path_to_output,
                verbose=verbose, write_file=True)
    geo_file = geo_inst.outfiles[0] + '.TAB'

    # Create instance with calibrated data

    cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst, 
                    verbose=verbose, write_file=True, 
                    local_path_to_output=global_path_to_output,
                    pnf_order=pnf_order, interact=False)
    cal_file = cal_inst.outfiles[0] + '.TAB'

    # Create instance with diffraction-limited profile and other
    # inputs needed for diffraction correction
    # (create_dlps returns None for missing ingress or egress)

    dlp_inst_ing, dlp_inst_egr = (
        rss.calibration.DiffractionLimitedProfile.create_dlps(
            rsr_inst, geo_inst, cal_inst, dr_km_desired,
            profile_range = profile_range,
            local_path_to_output=global_path_to_output,
            write_file=True, verbose=True))
    if direction == 'I':
        dlp_inst = dlp_inst_ing
    elif direction == 'E':
        dlp_inst = dlp_inst_egr
        
    dlp_file = dlp_inst.outfiles[0]+'.TAB'
    
    # create version of dlp instance required for diffraction reconstruction
    # This version reads the just-created GEO,CAL,DLP files to confirm they are valid

    if PerformInversion:
        if use_CSVs:
            data = rss.ExtractCSVData(geo_file, cal_file, dlp_file)
            #print('demo_e2e_event:',data.history)
        else:
            data = data_from_inst(geo_inst,cal_inst,dlp_inst,geo_file,cal_file_dlp_file)

    # create instance of retrieved optical depth profile

        tstart = time.time()
        tau_inst = rss.DiffractionCorrection(
           data, res_km, rng=inversion_range, resolution_factor=res_factor,
           psitype=psitype, wtype=wtype, verbose=verbose)
#        print(dir(tau_inst)) # need to get upated names of contents such as input_res

        if compute_tau_threshold_:
            tau_inst.tau_threshold_vals = compute_tau_threshold(cal_file,tau_inst)

        if compute_dr_pole:           
            rho_corr_pole_km_vals,junk,junk=\
                radius_correction_pole_from_tau_inst(dlp_file,tau_inst,
                        fit_number=fit_number,verbose=verbose)
            tau_inst.rho_corr_pole_km_vals = rho_corr_pole_km_vals
            
        if compute_dr_trajectory:
            rho_corr_timing_km_vals,junk,junk = \
                    radius_correction_trajectory_from_tau_inst(dlp_file,tau_inst,
                                    fit_number=fit_number,verbose=verbose)
            tau_inst.rho_corr_timing_km_vals =rho_corr_timing_km_vals

        tend = time.time()
        # print('TP1:tau_inst.tau_threshold_vals[0:10]',tau_inst.tau_threshold_vals[0:10])
        # print('TP1:tau_inst.rho_corr_pole_km_vals[0:10]',tau_inst.rho_corr_pole_km_vals[0:10])
        # print('TP1:tau_inst.rho_corr_timing_km_vals[0:10]',tau_inst.rho_corr_timing_km_vals[0:10])

        tau_file = write_tau_files_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
               res_km,inversion_range,res_factor,psitype,wtype,program,
                local_path_to_output=local_path_to_output)

        # plot results and save figfile
        figfile = None
        results = {'figfile':figfile,'geo_file':geo_file,'cal_file':cal_file,'dlp_file':dlp_file,
                   'tau_file':tau_file,'CORSS_8001_filepath':CORSS_8001_filepath,'res_km':res_km,
                   'res_factor':res_factor,'inversion_range':inversion_range,'psitype':psitype,
                   'wtype':wtype,'Rev':Rev,'direction':direction,'DSN':DSN[1:],'band':band,
                   'tau_threshold':tau_inst.tau_threshold_vals,'kernels':kernels}
        if plot:
            figfile = plot_comparisons(program,results,plot_phase=plot_phase,y_offset=y_offset,
                                       ylim_tau=ylim_tau,include_CORSS_8001=include_CORSS_8001,
                                       adjust_phase=adjust_phase,dy_offset=dy_offset,show=show,
                                      path_to_figs=global_path_to_demo_figs,save_figfile=save_figfile)
            results['figfile'] = figfile

        return(results)
    else:        
        data = rss.ExtractCSVData(geo_file, cal_file, dlp_file) # just for debugging
        return data,geo_file,cal_file,dlp_file,geo_inst,cal_inst,dlp_inst

def demo_e2e_event_loop(Rev='Rev007',direction='E',DSN='X43',
                        dr_km_desired=0.25,res_km=1.0,resolution_factor=0.75,
    profile_range=[85000., 90000.],inversion_range=[87475,87560],psitype='fresnel',
                        rsr_file = None,
                        _16kHz=None,
                decimate_16khz_to_8khz=False,
                decimate_16khz_to_4khz=False,
                decimate_16khz_to_2khz=False,
                decimate_16khz_to_1khz=False,
                decimate_50khz_to_25khz=False,
                decimate_50khz_to_10khz=False,
                decimate_50khz_to_5khz=False,
                decimate_50khz_to_1khz=False,
    wtype='kbmd20',include_CORSS_8001=True,plot=True,plot_phase=False,local_path_to_data = '../data/',
    local_path_to_tables='../tables/',save_figfile=True,
    adjust_phase=True,y_offset=0,dy_offset=0,local_path_to_kernels = '../kernels/',title=None,               
    verbose=True,silent=False,show=True,program='demo_e2e_event',use_CSVs=True,PerformInversion=True,
    compute_tau_threshold_ = True, compute_dr_pole=True, compute_dr_trajectory=True,fit_number=1):

#    print('ON ENTRY: psitype=',psitype)
    if _16kHz == None:
        if dr_km_desired <= 0.15:
            high_resolution = True # set to True for reconstructed resolution below 500 meters
        else:
            high_resolution = False
        _16kHz=high_resolution
        
    assert dr_km_desired<=min_required_dr_km_desired(res_km,resolution_factor),\
    'dr_km_desired '+str(dr_km_desired)+' is greater than '+str(min_required_dr_km_desired(res_km,resolution_factor))

#    verbose = True # set to False to suppress progress reports
#    silent = False # set to True to get confirmation that files exist, else False
#    show = True # set to False to suppress showing of plot at end (figfile still produced)
    
    planet = 'Saturn' 
    spacecraft = 'Cassini'
    pnf_order = 3 # polynomial order for normalizing power vs time
    CORSS_8001_filepath =get_CORSS_8001_XKa_TABfiles(CORSS_8001_all_filepaths=None,Rev=Rev,\
                            DSN=DSN,direction=direction,
                            local_path_to_data=local_path_to_data,silent=False)
    #print(CORSS_8001_filepath)

    if type(psitype) == type('fresnel'):
        psitypes = [psitype] # make it a list to enable loop below
    else:
        psitypes = psitype # a list already
    npsitypes = len(psitypes)
    
    if (type(res_km) == type(1.0)) or (type(res_km) == type(1)):
        res_kms = [res_km] # make it a list to enable loop below
    else:
        res_kms = res_km # a list already
    nres_kms = len(res_kms)

    if (type(resolution_factor) == type(0.75)) or (type(resolution_factor) == type(1)):
        resolution_factors=[resolution_factor]
    else:
        resolution_factors = resolution_factor
    nresolution_factors = len(resolution_factors)

    if type(wtype) == type('kmbd20'):
        wtypes=[wtype]
    else:
        wtypes = wtype
    nwtypes = len(wtypes)

# if multiple psitypes, make sure others are same length or scalars
    if npsitypes >1:
        if nresolution_factors == 1:
            resolution_factors = []
            for i in range(npsitypes):
                resolution_factors.append(resolution_factor)
        elif nresolution_factors != npsitypes:
            print("TROUBLE! incompatible numbers of psitypes and resolution_factors")
            return None
        if nres_kms ==1:
            res_kms = []
            for i in range(npsitypes):
                res_kms.append(res_km)
        elif nres_kms != npsitypes:
            print("TROUBLE! incompatible numbers of psitypes and res_kms")
            return
        if nwtypes ==1:
            wtypes = []
            for i in range(npsitypes):
                wtypes.append(wtype)
        elif nwtype != npsitypes:
            print("TROUBLE! incompatible numbers of psitypes and wtypes")
            return
# if multiple resolution_factors, make sure others are same length or scalars
    if nresolution_factors >1:
        if npsitypes == 1:
            psitypes = []
            for i in range(nresolution_factors):
                psitypes.append(psitype)
        elif nresolution_factors != npsitypes:
            print("TROUBLE! incompatible numbers of resolution_factors and psitypes")
            return None
        if nres_kms ==1:
            res_kms = []
            for i in range(nresolution_factors):
                res_kms.append(res_km)
        elif nres_kms != nresolution_factors:
            print("TROUBLE! incompatible numbers of resolution_factors and res_kms")
            return
        if nwtypes ==1:
            wtypes = []
            for i in range(nresolution_factors):
                wtypes.append(wtype)
        elif nwtypes != n_resfactors:
            print("TROUBLE! incompatible numbers of resolution_factors and wtypes")
            return
    # if multiple res_kms, make sure others are same length or scalars
    if nres_kms >1:
        if npsitypes == 1:
            psitypes = []
            for i in range(nres_kms):
                psitypes.append(psitype)
        elif nres_kms != npsitypes:
            print("TROUBLE! incompatible numbers of psitypes and res_kms")
            return None
        if nresolution_factors ==1:
            resolution_factors = []
            for i in range(nres_kms):
                resolution_factors.append(resolution_factor)
        elif nres_kms != nresolution_factors:
            print("TROUBLE! incompatible numbers of resolution_factors and res_kms")
            return
        if nwtypes ==1:
            wtypes = []
            for i in range(nres_kms):
                wtypes.append(wtype)
        elif nres_kms != nwtypes:
            print("TROUBLE! incompatible numbers of resolution_factors and wtypes")
            return
# if multiple wtypes, make sure others are same length or scalars
    if nwtypes >1:
        if npsitypes == 1:
            psitypes = []
            for i in range(nwtypes):
                psitypes.append(psitype)
        elif nwtypes != npsitypes:
            print("TROUBLE! incompatible numbers of psitypes and wtypes")
            return None
        if nresolution_factors ==1:
            resolution_factors = []
            for i in range(nwtypes):
                resolution_factors.append(resolution_factor)
        elif nwtypes != nresolution_factors:
            print("TROUBLE! incompatible numbers of resolution_factors and wtypes")
            return
        if nres_kms ==1:
            res_kms = []
            for i in range(nwtypes):
                res_kms.append(res_km)
        elif nres_kms != nwtypes:
            print("TROUBLE! incompatible numbers of res_km and wtypes")
            return

    # get required kernel files if not already fetched

    kernels = get_kernels_from_web(kernels_for_demo(),\
                local_path_to_kernels=local_path_to_kernels,silent=silent)
    kernel_naif = pick_SC_kernel(CORSS_8001_filepath)
    kernelSC = get_kernels_from_web(kernel_naif,local_path_to_kernels=local_path_to_kernels,silent=True)[0]
    kernels.append(kernelSC)

    # get required data files for all available DSNs if not already fetched
    if rsr_file == None:
        rsr_file = local_path_to_data + \
           pick_RSR_file(_16kHz=_16kHz,Rev=Rev,direction=direction,DSN=DSN,download=True,silent=True)

    # set up specifications for creating GEO/CAL/DLP files

    # Create instance with rsr file contents

    rsr_inst = rss.rsr_reader.RSRReader(rsr_file,verbose=True,
                                        decimate_16khz_to_8khz=decimate_16khz_to_8khz,
                                        decimate_16khz_to_4khz=decimate_16khz_to_4khz,
                                        decimate_16khz_to_2khz=decimate_16khz_to_2khz,
                                        decimate_16khz_to_1khz=decimate_16khz_to_1khz,
                                        decimate_50khz_to_25khz=decimate_50khz_to_25khz,
                                        decimate_50khz_to_10khz=decimate_50khz_to_10khz,
                                        decimate_50khz_to_5khz=decimate_50khz_to_5khz,
                                        decimate_50khz_to_1khz=decimate_50khz_to_1khz
                                       )

    # Create instance with geometry parameters

    geo_inst = rss.occgeo.Geometry(rsr_inst, planet, spacecraft,
                kernels,local_path_to_tables=global_path_to_tables,
                local_path_to_output=global_path_to_output,
                verbose=verbose, write_file=True)
    geo_file = geo_inst.outfiles[0] + '.TAB'

    # Create instance with calibrated data

    cal_inst = rss.calibration.Calibration(rsr_inst, geo_inst, 
                    verbose=verbose, write_file=True, 
                    local_path_to_output=global_path_to_output,
                    pnf_order=pnf_order, interact=False)
    cal_file = cal_inst.outfiles[0] + '.TAB'

    # Create instance with diffraction-limited profile and other
    # inputs needed for diffraction correction
    # (create_dlps returns None for missing ingress or egress)

    dlp_inst_ing, dlp_inst_egr = (
        rss.calibration.DiffractionLimitedProfile.create_dlps(
            rsr_inst, geo_inst, cal_inst, dr_km_desired,
            profile_range = profile_range,
            local_path_to_output=global_path_to_output,
            write_file=True, verbose=True))
    if direction == 'I':
        dlp_inst = dlp_inst_ing
    elif direction == 'E':
        dlp_inst = dlp_inst_egr
        
    dlp_file = dlp_inst.outfiles[0]+'.TAB'
    
    if PerformInversion:
    # create version of dlp instance required for diffraction reconstruction
    # This version reads the just-created GEO,CAL,DLP files to confirm they are valid

        if use_CSVs:
            data = rss.ExtractCSVData(geo_file, cal_file, dlp_file)
        else:
            data = data_from_inst(geo_inst,cal_inst,dlp_inst,geo_file,cal_file,dlp_file)
    
        results_list = []
        for psitype,resolution_factor,res_km,wtype in zip(psitypes,resolution_factors,res_kms,wtypes):
        
        # create instance of retrieved optical depth profile
            tstart = time.time()
            tau_inst = rss.DiffractionCorrection(
               data, res_km, rng=inversion_range, resolution_factor=resolution_factor,
               psitype=psitype, wtype=wtype, verbose=verbose)
            tend = time.time()

            if compute_tau_threshold_:
                tau_inst.tau_threshold_vals = compute_tau_threshold(cal_file,tau_inst)
    
            if compute_dr_pole:
                tau_inst.rho_corr_pole_km_vals = radius_correction_pole_from_tau_inst(dlp_file,tau_inst,fit_number=fit_number,verbose=verbose)
    
            if compute_dr_trajectory:
                tau_inst.rho_corr_timing_km_vals =radius_correction_trajectory_from_tau_inst(dlp_file,tau_inst,
                                                                                             fit_number=fit_number,verbose=verbose)

            # print('TP1:tau_inst.tau_threshold_vals[0:10]',tau_inst.tau_threshold_vals[0:10])
            # print('TP1:tau_inst.rho_corr_pole_km_vals[0:10]',tau_inst.rho_corr_pole_km_vals[0:10])
            # print('TP1:tau_inst.rho_corr_timing_km_vals[0:10]',tau_inst.rho_corr_timing_km_vals[0:10])
            
            tau_file = write_tau_files(tau_inst,geo_inst,cal_inst,dlp_inst,dlp_file,tstart,tend,
                           res_km,inversion_range,resolution_factor,psitype,wtype,program,
                            local_path_to_output = global_path_to_output)
        
         # plot results and save figfile
            figfile = None
            results = {'figfile':figfile,'rsr_file':rsr_file,'geo_file':geo_file,'cal_file':cal_file,'dlp_file':dlp_file,
                       'tau_file':tau_file,'CORSS_8001_filepath':CORSS_8001_filepath,'res_km':res_km,
                       'resolution_factor':resolution_factor,'inversion_range':inversion_range,'psitype':psitype,
                       'wtype':wtype,'Rev':Rev,'direction':direction,'DSN':DSN}
            if plot:
                figfile = plot_comparisons(program,results,plot_phase=plot_phase,y_offset=y_offset,
                          adjust_phase=adjust_phase,dy_offset=dy_offset,show=show,title=title,
                          save_figfile=save_figfile)
                results['figfile'] = figfile
                
            results_list.append(results)
        return(results_list)
    else:
        data = rss.ExtractCSVData(geo_file, cal_file, dlp_file) # just for debugging
        return data,geo_file,cal_file,dlp_file,geo_inst,cal_inst,dlp_inst

def demo_Rev007_X43E_Maxwell(show=False,_16kHz=False,
    program='demo_Rev007_X43E_Maxwell',
    Rev = 'Rev007',
    direction = 'E',
    DSN='X43',
    dr_km_desired=0.25,
    res_km=1.0,
    res_factor=0.75,
    profile_range=[85000., 90000.],
    inversion_range=[87475,87560],
    psitype='fresnel',
    wtype='kbmd20',
    include_CORSS_8001=True,
    plot=True,
    plot_phase=False,
    ylim_tau = (None,None),
    verbose=False,
    use_CSVs = True,
    silent=True):
    '''
    Demonstration of end-to-end processing for Maxwell ringlet Rev007_X43E
    '''
    
    results = demo_e2e_event(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,include_CORSS_8001=True,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,use_CSVs=use_CSVs,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results

def demo_Rev007_X43E_Maxwell_loop(show=False,_16kHz=False,
    decimate_16khz_to_1khz=True,
    decimate_16khz_to_2khz=False,
    decimate_16khz_to_4khz=False,
    decimate_16khz_to_8khz=False,
    program='demo_Rev007_X43E_Maxwell_loop',
    Rev = 'Rev007',
    direction = 'E',
    DSN='X43',
    dr_km_desired=0.25,
    res_km=[1.0,0.5],
    res_factor=0.75,
    profile_range=[85000., 90000.],
    inversion_range=[87475,87560],
    psitype='fresnel',
    wtype='kbmd20',
    include_CORSS_8001=True,
    plot=True,
    plot_phase=False,
    ylim_tau = (None,None),
    verbose=False,
    silent=True):
    '''
    Demonstration of end-to-end processing for Maxwell ringlet Rev007_X43E for multiple resolutions
    '''
    
    results = demo_e2e_event_loop(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,include_CORSS_8001=True,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,
        decimate_16khz_to_1khz=decimate_16khz_to_1khz,decimate_16khz_to_2khz=decimate_16khz_to_2khz,
        decimate_16khz_to_4khz=decimate_16khz_to_4khz,decimate_16khz_to_8khz=decimate_16khz_to_8khz,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results

def demo_Rev007_X43I_Maxwell(show=False,_16kHz=False,
    program='demo_Rev007_X43I_Maxwell',
    Rev = 'Rev007',
    direction = 'I',
    DSN='X43',
    dr_km_desired=0.25,
    res_km=1.0,
    res_factor=0.75,
    profile_range=[85000., 90000.],
    inversion_range=[87475,87560],
    psitype='fresnel',
    wtype='kbmd20',
    include_CORSS_8001=True,
    plot=True,
    plot_phase=False,
    ylim_tau = (None,None),
    verbose=False,
    use_CSVs=True,
    silent=True):
    '''
    Demonstration of end-to-end processing for Maxwell ringlet Rev007_X43I
    '''
    results = demo_e2e_event(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,include_CORSS_8001=True,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,use_CSVs=use_CSVs,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results
    
def demo_Fring(Rev = 'Rev054',direction = 'I',DSN='K55',show=False,dr_km_desired=0.075,psitype='fresnel',
               res_km=1.0,res_factor=0.75,profile_range=[138000.,142000.],inversion_range=[140040,140100.],
               wtype='kbmd20',include_CORSS_8001=True,silent=False,program='demo_Fring',
               ring_frame = 'IAU_SATURN',ylim_tau=(None,None),
               _16kHz=True,plot=True,plot_phase=False,adjust_phase=False,verbose=True):
    """
    Documentation for demo_Fring
    """
    results = demo_e2e_event(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,ring_frame=ring_frame,
        include_CORSS_8001=include_CORSS_8001,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,adjust_phase=adjust_phase,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results

def demo_StrangeRinglet(Rev = 'Rev067',direction = 'E',DSN='X14',show=False,dr_km_desired=0.075,psitype='newton',
               res_km=1.0,res_factor=0.75,profile_range=[108000.,128000.],inversion_range=[117800.,118000.],
               wtype='kbmd20',ylim_tau=(None,None),include_CORSS_8001=True,silent=False,program='demo_strange',
               ring_frame = 'IAU_SATURN',
               _16kHz=True,plot=True,plot_phase=False,adjust_phase=False,verbose=True):
    """
    Documentation for demo_StrangeRinglet
    """
    results = demo_e2e_event(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,ring_frame=ring_frame,
        include_CORSS_8001=include_CORSS_8001,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,adjust_phase=adjust_phase,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results

def demo_Rev133I_X25_Cripples_loop(show=False,
    _16kHz=True,
    decimate_16khz_to_1khz=True,
    decimate_16khz_to_2khz=False,
    decimate_16khz_to_4khz=False,
    decimate_16khz_to_8khz=False,
    program='demo_Rev133_X43E_Cripples_loop',
    Rev = 'Rev133',
    direction = 'I',
    DSN='X25',
    dr_km_desired = 0.025,
    res_km=[1.0,0.5],
    res_factor=0.75,
    profile_range=[65000., 85000.],
    inversion_range=[74400,77800],
    psitype='Newton',
    wtype='kbmd20',
    include_CORSS_8001=False,
    plot=True,
    plot_phase=False,
    ylim_tau = (None,None),
    verbose=False,
    silent=True):
    '''
    Demonstration of end-to-end processing for inner C ring ripple region Rev133_X25I for multiple resolutions
    '''
    
    results = demo_e2e_event_loop(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,include_CORSS_8001=True,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,
        decimate_16khz_to_1khz=decimate_16khz_to_1khz,decimate_16khz_to_2khz=decimate_16khz_to_2khz,
        decimate_16khz_to_4khz=decimate_16khz_to_4khz,decimate_16khz_to_8khz=decimate_16khz_to_8khz,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results
def demo_Rev125I_X34_Cripples_loop(show=False,
    _16kHz=True,
    decimate_16khz_to_1khz=True,
    decimate_16khz_to_2khz=False,
    decimate_16khz_to_4khz=False,
    decimate_16khz_to_8khz=False,
    program='demo_Rev125I_X34_Cripples_loop',
    Rev = 'Rev125',
    direction = 'I',
    DSN='X34',
    dr_km_desired = 0.025,
    res_km=[1.0,0.5],
    res_factor=0.75,
    profile_range=[65000., 85000.],
    inversion_range=[75100,75675],
    psitype='Newton',
    wtype='kbmd20',
    include_CORSS_8001=False,
    plot=True,
    plot_phase=False,
    ylim_tau = (None,None),
    verbose=False,
    silent=True):
    '''
    Demonstration of end-to-end processing for inner C ring ripple region RSS_125I_X34 for multiple resolutions
    '''
    
    results = demo_e2e_event_loop(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,include_CORSS_8001=True,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,
        decimate_16khz_to_1khz=decimate_16khz_to_1khz,decimate_16khz_to_2khz=decimate_16khz_to_2khz,
        decimate_16khz_to_4khz=decimate_16khz_to_4khz,decimate_16khz_to_8khz=decimate_16khz_to_8khz,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results

def demo_Rev007E_X43_W7493_loop(show=False,_16kHz=False,
    decimate_16khz_to_1khz=True,
    decimate_16khz_to_2khz=False,
    decimate_16khz_to_4khz=False,
    decimate_16khz_to_8khz=False,
    program='demo_Rev007_X43E_Maxwell_loop',
    Rev = 'Rev007',
    direction = 'E',
    DSN='X43',
    dr_km_desired=0.025,
    res_km=[1.0,0.5],
    res_factor=0.75,
    profile_range=[71000., 79000.],
    inversion_range=[74925,74945],
    psitype='fresnel',
    wtype='kbmd20',
    include_CORSS_8001=True,
    plot=True,
    plot_phase=False,
    ylim_tau = (None,None),
    verbose=False,
    silent=True):
    '''
    Demonstration of end-to-end processing for W74.93 wave Rev007_X43E for multiple resolutions
    '''
    
    results = demo_e2e_event_loop(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,include_CORSS_8001=True,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,
        decimate_16khz_to_1khz=decimate_16khz_to_1khz,decimate_16khz_to_2khz=decimate_16khz_to_2khz,
        decimate_16khz_to_4khz=decimate_16khz_to_4khz,decimate_16khz_to_8khz=decimate_16khz_to_8khz,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results

def demo_Rev133I_X25_W7494_loop(show=False,_16kHz=False,
    decimate_16khz_to_1khz=True,
    decimate_16khz_to_2khz=False,
    decimate_16khz_to_4khz=False,
    decimate_16khz_to_8khz=False,
    program='demo_Rev133I_X25_W7494_loop',
    Rev = 'Rev133',
    direction = 'I',
    DSN='X25',
    dr_km_desired=0.025,
    res_km=[1.0,0.5],
    res_factor=0.75,
    profile_range=[65000., 85000.],
    inversion_range=[74930,74950],
    psitype='Newton',
    wtype='kbmd20',
    include_CORSS_8001=False,
    plot=True,
    plot_phase=False,
    PerformInversion=True,
    ylim_tau = (None,None),
    verbose=False,
    rsr_file = None,
    CORSS_8001_filepath=None,
    silent=True):
    '''
    Demonstration of end-to-end processing for W74.93 wave Rev133I_X43 for multiple resolutions
    '''
    print("TP2: demo_Rev133I_X25_W7494_loop _16kHz=",_16kHz)
    
    results = demo_e2e_event_loop(Rev=Rev,direction=direction,DSN=DSN,
        dr_km_desired=dr_km_desired,res_km=res_km,res_factor=res_factor,
        profile_range=profile_range,inversion_range=inversion_range,psitype=psitype,
        wtype=wtype,include_CORSS_8001=True,plot=plot,plot_phase=plot_phase, 
        ylim_tau=ylim_tau,PerformInversion=PerformInversion,
        verbose=verbose,silent=silent,show=show,_16kHz=_16kHz,
        decimate_16khz_to_1khz=decimate_16khz_to_1khz,decimate_16khz_to_2khz=decimate_16khz_to_2khz,
        decimate_16khz_to_4khz=decimate_16khz_to_4khz,decimate_16khz_to_8khz=decimate_16khz_to_8khz,
        rsr_file = rsr_file, CORSS_8001_filepath=CORSS_8001_filepath,
        local_path_to_kernels=global_path_to_kernels,local_path_to_data=global_path_to_data,
        local_path_to_tables=global_path_to_tables,local_path_to_output=global_path_to_output)
    return results

def differential_opacity_plot(rev ='007',ID='Rev007E_rgf',direc ='E',
    profile = ['A inner','B4','A outer','Cas Div Ramp','B2','Maxwell','B1(B-flat)','C mid','C ramp'],
    dcol = ['blue','red','cyan','0.5','orange','0.0','crimson','lime','olive'],Lwhich=None,
    rho_min_vals = np.array([124000,112000,134000,121000,101500,87475,94500,82250,90750]),
    rho_max_vals = np.array([128000,115000,136000,122000,103000,87560,95250,83750,91800]),
    a_max_vals = [3e9,1e10],a_min_vals = np.logspace(6,9,301),
    q_vals = [2.8,3.0,3.2,3.4],m_vals=[1.78+0j],ylim=(-30,50),
    binsize=40.,edge_adjust=0,
    tau_min = 0.01,tau_max=4.0,Nmax=5000,TAUmatch='TAU_01KM.TAB'):
    '''
    Plot differential opacity models and observations
    '''
    if Lwhich == None:
        Lwhich = range(len(profile))
    if Lwhich != range(len(profile)):
        profile = [profile[i] for i in range(len(profile)) if i in Lwhich]
        dcol = [dcol[i] for i in range(len(dcol)) if i in Lwhich]
        rho_min_vals = [rho_min_vals[i] for i in range(len(rho_min_vals)) if i in Lwhich]
        rho_max_vals = [rho_max_vals[i] for i in range(len(rho_max_vals)) if i in Lwhich]
    #print('debug:profile',profile)    
    wildcard = global_path_to_data + 'CORSS_8001/data/Rev'+rev+'/Rev*'+direc+'/*/*'+TAUmatch
    
    # END OF USER INPUT
    
    q_vals,m_vals,tau_ratio_KX,tau_ratio_XS,markers,lines,colors,labels = \
        particle_size_models(q_vals=q_vals,m_vals=m_vals)
    
    paths = glob.glob(wildcard,flags=glob.BRACE)
    files = [os.path.basename(pathi) for pathi in paths]
    
    # storage lists
    dc = []
    nfiles = len(files)
    
    nprofiles = len(profile)
    tau_k = np.full((nprofiles,Nmax),np.nan)
    tau_x = np.full((nprofiles,Nmax),np.nan)
    tau_s = np.full((nprofiles,Nmax),np.nan)
    
    # read in data and bin
    for ifile in range(nfiles):
    
        rho,tau = np.loadtxt(paths[ifile],delimiter=',',usecols=(0,6)).T
    
        # iterate over all profiles
        jprofile = 0
        for rho_min,rho_max,p in zip(rho_min_vals,rho_max_vals,profile):
    
            tau_vals = []
    
            # bin up data in binsize km bins, adjusting edges by edge_adjust to avoid zero contributions near edges
            for rhoi in np.arange(rho_min+edge_adjust,rho_max-edge_adjust,binsize):
    
                # bin within limits and within binsize km of bin center
                bins = [(rho>rhoi-binsize/2)&(rho<rhoi+binsize/2)&(rho>rho_min)&(rho<rho_max)][0]
    
                # store
                if len(bins)>0:
                    # idiom to set taui to np.nan if tau[bins] is entirely nan
                    taui = np.nan if np.all(tau[bins]!=tau[bins]) else np.nanmean(tau[bins])
                    if taui != np.nan:
                        tau_vals.append(taui) #tau_vals += [taui]
            # store in band-specific list
            ntau = len(tau_vals)
            if ntau >0:
                if '_K' in files[ifile]:
                    istart = find_first_nan_index(tau_k[jprofile,:])
        #            print('istart=',istart)
                    tau_k[jprofile,istart:istart+ntau] = tau_vals
                elif '_X' in files[ifile]:
                    istart = find_first_nan_index(tau_x[jprofile,:])
                    tau_x[jprofile,istart:istart+ntau] = tau_vals
                elif '_S' in files[ifile]:
                    istart = find_first_nan_index(tau_s[jprofile,:])
                    tau_s[jprofile,istart:istart+ntau] = tau_vals
            jprofile += 1
    
    # d tau / tau
    dt_kx = 100.*(tau_k-tau_x)/tau_x
    dt_xs = 100.*(tau_x-tau_s)/tau_x
    
    '''
        ~~ Plot results! ~~
    '''
    plt.plot([-100,100],[0,0],dashes=[12,4],color='0.5')
    plt.plot([0,0],[-100,100],dashes=[12,4],color='0.5')
    out = open('dtau_'+ID+'_rgf.csv','w')
    headers = ['Feature','rho_min (km)','rho_max (km)','tau_kx','tau_xs','sigma_tau_kx','sigma_tau_xs','sigma_total']
    out.write(",".join(str(item) for item in headers)+'\n')
    # loop over features
   # print('rho_min_vals',rho_min_vals)
    
    for jprofile in np.argsort(rho_min_vals):
    #    print('profile',profile[jprofile])
        #print('profile,jprofile,profile[jprofile]',profile,jprofile,profile[jprofile])
        this_label = profile[jprofile]
        #print('this_label',this_label)
        L=np.where((tau_k[jprofile,:]>tau_min)&(tau_k[jprofile,:]<tau_max)&(tau_x[jprofile,:]>tau_min)&
                   (tau_x[jprofile,:]<tau_max)&(tau_s[jprofile,:]>tau_min)&(tau_s[jprofile,:]<tau_max))[0]
    #    print('len(L)',len(L))
        if len(L) >0 :
            plt.plot(dt_kx[jprofile,L],dt_xs[jprofile,L],'.',color=dcol[jprofile],zorder=4,label=this_label)
    #        print(profile[jprofile],dt_kx[jprofile,L],dt_xs[jprofile,L])
        # highlight data with a circle
            kxm = np.nanmedian(dt_kx[jprofile,L])
            xsm = np.nanmedian(dt_xs[jprofile,L])
            kxs = np.sqrt(np.nanmedian(np.square(dt_kx[jprofile,L]-kxm)))
            xss = np.sqrt(np.nanmedian(np.square(dt_xs[jprofile,L]-xsm)))
            rad = np.sqrt(kxs**2.+xss**2.)/np.sqrt(2)
            row = [profile[jprofile],rho_min_vals[jprofile],rho_max_vals[jprofile],round(kxm,3),round(xsm,3),round(kxs,3),round(xss,3),round(rad,3)]
            out.write(",".join(str(item) for item in row)+'\n')
            ekg = {'xy':[kxm,xsm],'width':5.*kxs,'height':5.*xss,'color':dcol[jprofile],'fill':False,'zorder':5}
            elps = pch.Ellipse(**ekg)
            plt.gcf().gca().add_artist(elps)
        else:
            print('Skipping',profile[jprofile])
    out.close()
    # text label offsets
    xoff = [-10,-10,-2,-15,3,3,3]
    yoff = [-4,-4,-4,-2,0,0,0]
    
    # iterate over all model power laws
    for m in m_vals:
        for i in range(len(labels)):
    
        # model grid
            if i == 0 or i == 5:
                pkwargs = {'color':colors[i],'ls':lines[i],'label':labels[i],'lw':1}
            else:
                pkwargs = {'color':colors[i],'ls':lines[i],'lw':1}
            plt.plot(100.*tau_ratio_KX[i],100.*tau_ratio_XS[i],**pkwargs)
        
            # q label
            if i < int(len(a_max_vals)*len(q_vals)/2):
                tx = 100.*tau_ratio_KX[i][0]+5
                ty = 100.*tau_ratio_XS[i][0]-2
                plt.text(tx,ty,'q='+str(q_vals[i]),clip_on=True)
        
            # a_min labels
            if i == np.min([len(labels)-1,3]):
                pkwargs = {'color':colors[i],'ls':lines[i],'label':labels[i],'lw':1}
                for j in range(0,len(tau_ratio_KX[i]),50):
                    if j < 100:
                        lt = str(int(a_min_vals[j]/1e6))+' mm'
                    elif j >= 100 and j < 200:
                        lt = str(int(a_min_vals[j]/1e7))+' cm'
                    elif j >= 200 and j < 300:
                        lt = str(int(a_min_vals[j]/1e8)*int(10))+' cm'
                    else:
                        lt = str(int(a_min_vals[j]/1e9))+' m'
                    tx = 100.*tau_ratio_KX[i][j]+xoff[int(j/50)]
                    ty = 100.*tau_ratio_XS[i][j]+yoff[int(j/50)]
                    plt.text(tx,ty,lt,zorder=5,clip_on=True)
        
            # a_min markers
            pkwargs = {'color':'none','ls':'none','marker':markers[i],'markeredgecolor':colors[i]}
            for j in range(0,len(tau_ratio_KX[i]),50):
                plt.plot(100.*tau_ratio_KX[i][j],100.*tau_ratio_XS[i][j],**pkwargs)
        
        plt.xlim(-30,100)
        plt.ylim(ylim)
        plt.xlabel(r'$\Delta\tau_{Ka-X} / \tau_X$ as %')
        plt.ylabel(r'$\Delta\tau_{X-S} / \tau_X$ as %')
        plt.legend(loc=4,frameon=False,numpoints=1)
        plt.title(r'Scattering Model with $\overline{m}=$'+str(m))#+' for '+profile)
        figfile = global_path_to_local_figs + 'dtau_dtau_'+str(m)+'_'+ID+'.png'
        plt.savefig(figfile)
        print(figfile)
        plt.show()

# integrand of Delta tau
def dtau_int(a,dQ,q=3.1):
    return a**2. * dQ * n(a,q)

def event_name_from_rev_info(rev_info): # ex:  'RSS_133E_X43' as used in C ring ripples paper  
    event_name = 'RSS_'+rev_info['rev_num']+rev_info['prof_dir'][1]+'_'+rev_info['band'][1]+rev_info['dsn']
    return event_name

def event_name_from_taufilepath(taufilepath):
    rev_info = get_rev_info_from_tau(taufilepath)
    return event_name_from_rev_info(rev_info)

    
def featurelist_to_dict(infile = global_path_to_local_tables + 'my_feature_list_short_Cring.csv'):
# NB this fails for full feature list with ID in first column such as 14a in venv2 kernel
    dict_features = np.genfromtxt(infile, delimiter= ',',dtype=None,names=True,encoding=None)
    dict_features = {name: dict_features[name].tolist() for name in dict_features.dtype.names}
    dict_features = {"Baillie" if k == '\ufeffBaillie' else k:v for k,v in dict_features.items()} #get rid of '\ueff' in first label
    dict_features['a_res'] = [float(item) for item in dict_features['a_res']]
    dict_features['dr_min'] = [float(item) for item in dict_features['dr_min']]
    dict_features['dr_max'] = [float(item) for item in dict_features['dr_max']]
    return dict_features

def find_first_nan_index(arr):
  """
  Finds the index of the first NaN element in a NumPy array.

  Args:
    arr: A NumPy array.

  Returns:
    The index of the first NaN element, or None if no NaN is found.
  """
  nan_mask = np.isnan(arr)
  nan_indices = np.where(nan_mask)[0]
  if nan_indices.size > 0:
    return nan_indices[0]
  else:
    return None

def find_index(my_list, my_string,verbose=False):
    '''
    Find index of my_string in my_list
    '''
    this_list = my_list#.tolist()
#    print(this_list)
    try:
        index_value = this_list.index(my_string)
        if verbose:
            print('find_index: index_value for ',my_string,'is',index_value)
        return index_value
    except ValueError:
        print('Unable to find index value for',my_string)
        return -1
        
def flatten_extend(matrix):
    '''
    Flattens a matrix into a 1D list
    '''
    flat_list = []
    for row in matrix:
        flat_list.extend(row)
    return flat_list
    
def format_time(name,rev,band,dsn,direc,dlp_res_used,psitype,res_km,tstart,tend,MINUTES=False):
    '''
    Construct formatted execution time
    '''
    dt = tend - tstart
    if MINUTES:
        dt /= 60.
        units = ' minutes'
    else:
        units = ' seconds'
    return name+ ' RSS_'+rev+direc+'_'+band+dsn+' res '+f'{res_km:0.3f} km '+f'{psitype:13}'+dlp_res_used+' ' + f'time {dt:0.2f}'+units

def get_all_RSR_files(local_path_to_tables=global_path_to_tables,force=False,download=True,silent=False):
    '''
    Grab all RSR files from a specified list
    Warning: this requires high internet bandwidth and several hours
    '''
    file_list = local_path_to_tables+'rsr_16kHz_files_before_USO_failure.txt'
    with open(file_list,'r') as f:
        RSR_files = f.read().splitlines()
    for rsr_file in RSR_files:
        get_RSRfiles_from_web(rsr_file,force=force,silent=silent)
    file_list = local_path_to_tables+'rsr_all_files_before_USO_failure.txt'
    with open(file_list,'r') as f:
        RSR_files = f.read().splitlines()
    for rsr_file in RSR_files:
        get_RSRfiles_from_web(rsr_file,force=force,silent=silent)

def get_all_SC_kernels(kernels_list = global_path_to_tables+'e2e_kernels.ker',
              local_path_to_kernels=global_path_to_kernels,silent=False):
    '''
    Grab all Cassini SC kernels from kernels_list
    '''
    with open(kernels_list,'r') as f:
        kernels_all  = f.read().splitlines() # strips \n at end of each line.
    for kernel in kernels_all:
        try: # skip over lines in kernels_all that aren't Cassini S/C ephemerides
            index= kernel.index('SCPSE')
            index = kernel.index('$A/') +len('$A/')
            kernel_found = kernel[index:index+len('050623R_SCPSE_05132_05150.bsp')]
            get_kernels_from_web('naif/CASSINI/kernels/spk/'+kernel_found,
                                 local_path_to_kernels=local_path_to_kernels,silent=silent)
        except: # this line does not contain a Cassini S/C ephemeris
            continue

def get_CORSS_8001_file(CORSS_8001_filepath,local_path_to_data=global_path_to_data,force=False,silent=False): 
    '''
    Grab Cassini PDS RSS files from web
    '''
    webURL = 'https://pds-rings.seti.org/holdings/volumes/CORSS_8xxx/'
    return get_files_from_web(CORSS_8001_filepath,local_path_to_data,webURL,force=force,silent=silent)[0]

def get_CORSS_8001_TABfiles(CORSS_8001_all_filepaths=None,Rev=None,\
                            DSN=None,direction=None,local_path_to_data=global_path_to_data,
                                local_path_to_tables = global_path_to_tables,silent=True):
    '''
    Retrieve requested existing PDS3 RSS *TAU*.TAB files
    '''
    if CORSS_8001_all_filepaths==None:
        CORSS_file_list = local_path_to_tables+'CORSS_8001_TAU_filepaths_all.txt' # contains complete list of TAU TAB files in CORSS_8001 on PDS
        with open(CORSS_file_list,"r") as f:
            CORSS_8001_all_filepaths  = f.read().splitlines() # strips \n at end of each line
    if DSN=='K47' and Rev == 'Rev064':
        print('approximate match found')
        DSN='X43' # K47 not on CORSS and we just want to get the correct S/C kernel for this event:
    if DSN=='K34' and Rev == 'Rev064':
        print('approximate match found')
        DSN='X43' # K34 not on CORSS and we just want to get the correct S/C kernel for this event:
    if DSN=='K47' and Rev == 'Rev081':
        print('approximate match found')
        DSN='X43' # K47 not on CORSS and we just want to get the correct S/C kernel for this event:
    # if DSN=='X63' and Rev == 'Rev125':
    #     print('approximate match found')
    #     DSN='X43' # K47 not on CORSS and we just want to get the correct S/C kernel for this event:

    if DSN != None and direction != None:
        DSNdirection = DSN+'_'+direction
    else:
        DSNdirection = None
    filepaths = []
    for filepath in CORSS_8001_all_filepaths:
        if '#' in filepath: # skip comment lines
            continue
#        print('Rev:',Rev,'DSNdirection',DSNdirection)
        if (Rev == None or Rev in filepath) and (DSNdirection == None or DSNdirection in filepath):
            get_CORSS_8001_file(filepath,local_path_to_data=local_path_to_data,silent=silent)
            if not silent:
                print("filepath:",filepath)
            filepaths.append(filepath)
    if len(filepaths) == 1:
        return filepaths[0]
    else:
        return filepaths
        

def get_CORSS_8001_XKa_TABfiles(CORSS_8001_all_filepaths=None,Rev=None,\
                            DSN=None,direction=None,local_path_to_data=global_path_to_data,
                                local_path_to_tables = global_path_to_tables,silent=True):
    if CORSS_8001_all_filepaths==None:
        CORSS_file_list = local_path_to_tables+'CORSS_8001_TAU_filepaths.txt' # contains complete list of TAU TAB files in CORSS_8001 on PDS
        with open(CORSS_file_list,"r") as f:
            CORSS_8001_all_filepaths  = f.read().splitlines() # strips \n at end of each line
    if DSN != None and direction != None:
        DSNdirection = DSN+'_'+direction
    else:
        DSNdirection = None
    filepaths = []
    for filepath in CORSS_8001_all_filepaths:
        if '#' in filepath: # skip comment lines
            continue
        if (Rev == None or Rev in filepath) and (DSNdirection == None or DSNdirection in filepath):
            get_CORSS_8001_file(filepath,local_path_to_data=local_path_to_data,silent=silent)
            filepaths.append(filepath)
    if len(filepaths) == 1:
        return filepaths[0]
    else:
        return filepaths

def get_CORSS_8001_TAUfile(rev_info,local_path_to_data = global_path_to_data):
    Rev = rev_info['rev_num']
    DSN = rev_info['band'][1]+rev_info['dsn'][4:]
    direction = rev_info['prof_dir'][1]
    taufilepath = get_CORSS_8001_TABfiles(Rev=Rev,DSN=DSN,direction=direction)
    if len(taufilepath) == 0:
        print('No CORSS_8001 file for ',Rev,DSN,direction)
        print('taufilepath:',taufilepath)
        return None
    else:
        return local_path_to_data+taufilepath

def get_dBHz(_16kHz=True,Rev='007',direction='I',DSN='X43',
             snr0_file='snr0_rsr_all_files_before_USO_failure.csv',
                  local_path_to_tables=global_path_to_tables,verbose=False):
    '''
    Get SNR in dBHz by matching rsr filename to tabulated results
    '''
    rsr_file = pick_RSR_file(_16kHz=True,Rev=Rev,direction=direction,DSN=DSN,
                  local_path_to_tables=local_path_to_tables,force=False,download=False,silent=True)
    f = open(local_path_to_tables + snr0_file, "r")
    contents = f.read().splitlines() # strips \n at end of each line
    f.close()
    dBHz = np.nan
    for line in contents:
        if rsr_file in line:
            dBHz = float(line[line.index(',')+1:])
            if verbose:
                print('Found snr0=',dBHz,'for',os.path.basename(rsr_file))
            break
    return dBHz 

def get_dBHz_from_rsr_file(rsr_file,
             snr0_file='snr0_rsr_all_files_before_USO_failure.csv',
                  local_path_to_tables=global_path_to_tables,verbose=False):
    '''
    Get SNR in dBHz by matching rsr filename to tabulated results
    '''
    
    f = open(local_path_to_tables + snr0_file, "r")
    contents = f.read().splitlines() # strips \n at end of each line
    f.close()
    dBHz = np.nan
    for line in contents:
        if rsr_file.upper() in line.upper():
            dBHz = float(line[line.index(',')+1:])
            if verbose:
                print('Found snr0=',dBHz,'for',os.path.basename(rsr_file))
            break
    return dBHz 
       
def get_files_from_web(files,local_path,webURL,force=False,silent=True):
    '''
    Utility routine to get requested files from web URL
    '''
    local_files = []
    if type(files) is str: # if this is a single explicit string, convert to a list
        files = [files]
    for file in files:
        localdirfile = local_path + file
        local_files.append(localdirfile)
        webdirfile = webURL + file
        if not force and os.path.isfile(localdirfile):
            if not silent:
                print(localdirfile + ' already exists')
            continue
        elif force and os.path.isfile(localdirfile):
            print(localdirfile + ' already exists -- possibly short -- and is being overwritten...')
        print("Current file",localdirfile," Fetching...")
        cmd = 'curl --create-dirs -o '+localdirfile+' "'+webdirfile+'"'
        response = os.system(cmd)
        if os.path.isfile(localdirfile):
            if not silent:
                print(localdirfile + ' now exists')
        else:
            print("TROUBLE! ",localdirfile,' does not exist')
    return local_files


def get_kernels_from_web(kernels,local_path_to_kernels=global_path_to_kernels,force=False,silent=False):
    '''
    Grab kernels from web.
    NO check to seee if local/ is in the name!
    '''
    webURL = "https://naif.jpl.nasa.gov/pub/"
    return get_files_from_web(kernels,local_path_to_kernels,webURL,force=force,silent=silent)
    
def get_loaded_kernels(basename=True):
    count = spice.ktotal('ALL')
    files = []
    for i in range(count):
    
        file, ftype, srcfil, handle =   spice.kdata( i, 'ALL')
        if basename:
            files.append(os.path.basename(file))
        else:
            files.append(file)
    return files
    
def get_processor_info():
    '''
    Get processor information for several possible system platforms
    '''
    if platform.system() == "Windows":
        return platform.processor()
    elif platform.system() == "Darwin":
        info = subprocess.check_output(['/usr/sbin/sysctl', "-n", "machdep.cpu.brand_string"]).strip()
        return info.decode()
    elif platform.system() == "Linux":
        command = "cat /proc/cpuinfo"
        return subprocess.check_output(command, shell=True).strip()
    return ""

def get_psitype(labelfile):
    with open(labelfile, 'r') as file:
        contents = file.read().splitlines()
    for line in contents: # look for psitype in history part of label file
        if 'psitype: ' in line:
            line = line.strip()
            line = line.replace('psitype: ','')
            psitype = line.replace(',','')
            return(psitype.lower())
    return('')

def get_rev_info_from_tau(taufile,rsr_file='"UNKNOWN"'):
    if rsr_file != '"UNKNOWN"':
        rsr_file = '"' + os.path.basename(rsr_file) + '"'
    basename = os.path.basename(taufile)
    band = '"'+basename[basename.index('RSS_')+13]+'"'
    year = basename[basename.index('RSS_')+4:basename.index('RSS_')+8]
    doy  = basename[basename.index('RSS_')+9:basename.index('RSS_')+12]
    dsn  = basename[basename.index('RSS_')+14:basename.index('RSS_')+16]
    occ_dir = '"BOTH"' # not correct, possibly
    fields = taufile.split('/')
    #for i,f in enumerate(fields):
    #    print(i,f)
    direc = fields[-2][-1]
    rev_num = fields[-4][3:]
    #print(direc,rev_num)
    assert (direc == 'I' or direc == 'E'),"Illegal direction "+direc
    
    if direc == 'E':
        prof_dir = '"EGRESS"'
    elif direc == 'I':
        prof_dir = '"INGRESS"'
    planetary_occ_flag = '"Y"'
    rev_info = {'rsr_file':rsr_file,'band':band,'year':year,'doy':doy,'dsn':dsn,
               'occ_dir':occ_dir,'prof_dir':prof_dir,'rev_num':rev_num,'planetary_occ_flag':planetary_occ_flag}
    #print('\n',rev_info)
    return rev_info
    
def get_RSRfiles_from_web(RSRfiles,local_path_to_data=global_path_to_data,force=False,silent=False):
    '''
    Grab RSR files from web
    '''
    webURL = web="https://atmos.nmsu.edu/pdsd/archive/data/"
    return get_files_from_web(RSRfiles,local_path_to_data,webURL,force=force,silent=silent)

def get_trajectory_correction_coefficients(rev_info,fit_number=1,verbose=False,local_path_to_local_tables=global_path_to_local_tables):

    assert (fit_number==1 or fit_number==7),'radius_correction_pole: illegal fit_number '+str(fit_number)
    
    if fit_number == 1:
        fit_output_file = local_path_to_local_tables + 'ringfit_v1.8.Sa025S-RF-V5574.out'
    elif fit_number == 7: # weighted fit with multiple entries, final one overwrites previous, which is correct
        fit_output_file = local_path_to_local_tables + 'ringfit_v1.8.Sa025S-RF-V5981.out'

    with open(fit_output_file,'r') as f:
        contents = f.readlines()

    rev_num = rev_info['rev_num']
    direc = rev_info['prof_dir'][1]
    band = rev_info['band'][1]
    dsn = rev_info['dsn']

    match = 'RSS_'+rev_num+direc+'_'+band+dsn

    dt = 0.
    alpha = 0.
    for line in contents:
        if match in line:
            if 'D'+dsn+'_dt' in line:
                if verbose:
                    print(line)
                    print(line.split()[4])
                dt = float(line.split()[4])
            if 'D'+dsn+'_dr/dR1000' in line:
                if verbose:
                    print(line)
                    print(line.split()[4])
                alpha = float(line.split()[4])
    r0 = 100000.        
    return dt,alpha,r0

def is_res_km_valid(dlp_file,res_factor,res_km):
    '''
    Raise exception if requested res_km too small for given DLP file
    '''
    min_valid_res_km_ = min_valid_res_km(dlp_file,res_factor)
    if res_km < min_valid_res_km_:
        raise Exception(f"Requested res_km {res_km} is less than minimum allowed value {min_valid_res_km_:0.3f} for {dlp_file}")
        
def kernel_is_loaded(kernel):
    kernels_loaded = get_loaded_kernels(basename=True)
    if os.path.basename(kernel) in kernels_loaded:
        return True
    else:
        return False
        
def kernels_for_demo(verbose=False): 
    ''' 
    Grab Rev007 and F ring and Strange ringlet kernels
    '''
    kernels = [
        'naif/CASSINI/kernels/spk/050606R_SCPSE_05114_05132.bsp', # Rev007 
        'naif/CASSINI/kernels/pck/cpck26Feb2009.tpc',
        'naif/CASSINI/kernels/lsk/naif0012.tls',
        'naif/generic_kernels/spk/planets/de430.bsp',
        'naif/generic_kernels/spk/stations/a_old_versions/earthstns_itrf93_050714.bsp',
        'naif/CASSINI/kernels/pck/earth_000101_180919_180629.bpc',
        'naif/CASSINI/kernels/fk/earth_topo_050714.tf',
        'naif/generic_kernels/spk/stations/prelim/dss_47_prelim_itrf93_240908.bsp',
        'naif/CASSINI/kernels/fk/stations/prelim/dss_47_prelim_itrf93_240908.tf',
        'local/Saturn_F_ring_frame.tf',
        'local/Saturn_StrangeRinglet_20250605.tf',
        'local/Saturn_ring_plane.tf'
    ]
    if verbose:
        print('kernels to be loaded:',kernels)
    return kernels
    
def load_kernels_for_taufile(taufile,silent=True,load_absent_only=True): 
    '''
    Load the necessary kernels to compute values contained in a taufile
    '''
    rev_info = get_rev_info_from_tau(taufile)
    
    Rev=rev_info['rev_num']
    DSN=rev_info['dsn']
    direction = rev_info['prof_dir'][1]
    CORSS_8001_filepath =get_CORSS_8001_TABfiles(CORSS_8001_all_filepaths=None,Rev=Rev,\
                    DSN=DSN,direction=direction,local_path_to_data=global_path_to_data,silent=silent)[0]
#    print('CORSS_8001_filepath',CORSS_8001_filepath)
    kernels = get_kernels_from_web(kernels_for_demo(),\
                local_path_to_kernels=global_path_to_kernels,silent=silent)
    kernel_naif = pick_SC_kernel(CORSS_8001_filepath)
    kernelSC = get_kernels_from_web(kernel_naif,
                                    local_path_to_kernels=global_path_to_kernels,silent=silent)[0]
    kernels.append(kernelSC)

    for kernel in kernels:
        if load_absent_only and kernel_is_loaded(kernel):
            continue
        else:
            spice.furnsh(kernel)

def merge_sorted_tabfiles(indir,inversion_range,psitype=None,tmp2output=True,
                          tmp2globaloutput=False,add_inversion_range=False,update_radius_corrections = False,
                          clean=False,silent=False):
    firstfile = glob.glob(indir+'R*.TAB.*')[0] # assumes ALL files in this directory are to be merged!
    outfile = firstfile[0:-5]
    if not silent:
        print('nominal outfile',outfile)
    files = glob.glob(indir+'R*.TAB.*')
    r0 = np.zeros(len(files)) # sort by radius of first entry
    for i,file in enumerate(files):
        # print(file)
        d = np.loadtxt(file,delimiter=',')
#        print(i,d[0,0])
        r0[i] = d[0,0]
    index=np.argsort(r0)
#    outfile = indir + 'testit.TAB'
    os.system('cp /dev/null '+outfile) # create zero-length file
#    print(index)
    with open(outfile, 'a') as destination_file:
        for j in index:
            with open(files[j], 'r') as source_file:
                source_content = source_file.read()
            destination_file.write(source_content)
    if not silent:
        print('wrote ',outfile)
    with open(outfile, 'r') as file:
        contents = file.read()
        file_records = contents.count('\n')
    if not silent:
        print("Total lines:", file_records)
    if update_radius_corrections:
        if not silent:
            print("About to update radius corrections...")
            load_kernels_for_taufile(outfile,silent=True,load_absent_only=True)
            update_taufile_radius_corrections(outfile,update_dr_pole=True,update_dr_trajectory=True,
                                     fit_number=1,NMAX_POLE = 500,NMAX_TRAJECTORY = 1000,
                                     write_only_if_zero=True,verbose=not silent,
                                     save_original=True)
    data = np.loadtxt(outfile,delimiter=',')
    ring_radius = data[:,0]
    ring_longitude = data[:,3]
    ring_azimuth = data[:,4]
    oet = data[:,9]
    ret = data[:,10]
    scet = data[:,11]
    elevation = data[:,12]

    basename = os.path.basename(outfile)
    year = basename[4:8]
    doy = basename[9:12]

# construct label file
    first_tab_file = (files[index[0]])
    if not silent:
        print('first_tab_file:',first_tab_file)
    first_label_file = first_tab_file.replace('TAB','LBL')
    if not silent:
        print('first_label_file:',first_label_file)
    with open(first_label_file,'r') as f:
        contents_first = f.readlines()
    if psitype == None:
        psitype = get_psitype(first_label_file)
    last_tab_file = files[index[-1]]
    last_label_file = last_tab_file.replace('TAB','LBL')
    if not silent:
        print(last_label_file)
    with open(last_label_file,'r') as f:
        contents_last = f.readlines()

    iline = 0
    label_out = np.copy(contents_first)
# all newly-formatted lines need a trailing space that is later stripped since incoming lines have newline
    for line_first,line_last in zip(contents_first,contents_last):
        line_out = label_out[iline]
        if line_first != line_last:
            if not silent:
                print(line_first)
                print(line_last)
                print('\n')
        if 'FILE_RECORDS' in line_out:
            label_out[iline] = 'FILE_RECORDS                         = '+str(file_records)+' '
        if 'ROWS' in line_out:
            label_out[iline] = '  ROWS                        = '+str(file_records)+' '
        if 'MINIMUM_RING_RADIUS' in line_out:
            label_out[iline] = 'MINIMUM_RING_RADIUS                  = '+f'{np.min(ring_radius):10.3f}   <km> '
        if 'MAXIMUM_RING_RADIUS' in line_out:
            label_out[iline] = 'MAXIMUM_RING_RADIUS                  = '+f'{np.max(ring_radius):10.3f}   <km> '
        if 'MINIMUM_SAMPLING_PARAMETER' in line_out:
            label_out[iline] = '  MINIMUM_SAMPLING_PARAMETER  = '+f'{np.min(ring_radius):10.3f} '
        if 'MAXIMUM_SAMPLING_PARAMETER' in line_out:
            label_out[iline] = '  MAXIMUM_SAMPLING_PARAMETER  = '+f'{np.max(ring_radius):10.3f} '        
        if 'MINIMUM_RING_LONGITUDE' in line_out:
            label_out[iline] = 'MINIMUM_RING_LONGITUDE               = '+f'{np.min(ring_longitude):8.4f}   <deg> '
        if 'MAXIMUM_RING_LONGITUDE' in line_out:
            label_out[iline] = 'MAXIMUM_RING_LONGITUDE               = '+f'{np.max(ring_longitude):8.4f}   <deg> '
        if 'MINIMUM_OBSERVED_RING_AZIMUTH' in line_out:
            label_out[iline] = 'MINIMUM_OBSERVED_RING_AZIMUTH        = '+f'{np.min(ring_azimuth):8.4f}   <deg> '
        if 'MAXIMUM_OBSERVED_RING_AZIMUTH' in line_out:
            label_out[iline] = 'MAXIMUM_OBSERVED_RING_AZIMUTH        = '+f'{np.max(ring_azimuth):8.4f}   <deg> '
        if '^SERIES' in line_out:
            label_out[iline] = '^SERIES                              = "' + os.path.basename(outfile) + '" '
        if 'PRODUCT_ID                           =' in line_out:
            label_out[iline] = 'PRODUCT_ID                           = "' + os.path.basename(outfile) + '" '
        if 'START_TIME                           =' in line_out:
            label_out[iline] = 'START_TIME                           = '+ spm2date(year,doy,oet[0]) +' ' 
        if 'STOP_TIME                            =' in line_out:
            label_out[iline] = 'STOP_TIME                            = '+ spm2date(year,doy,oet[-1]) +' ' 
        if 'RING_EVENT_START_TIME                =' in line_out:
            label_out[iline] = 'RING_EVENT_START_TIME                = '+ spm2date(year,doy,ret[0]) +' ' 
        if 'RING_EVENT_STOP_TIME                 =' in line_out:
            label_out[iline] = 'RING_EVENT_STOP_TIME                 = '+ spm2date(year,doy,ret[-1]) +' ' 
        if 'SPACECRAFT_EVENT_START_TIME          =' in line_out:
            label_out[iline] = 'SPACECRAFT_EVENT_START_TIME          = '+ spm2date(year,doy,scet[0]) +' ' 
        if 'SPACECRAFT_EVENT_STOP_TIME           =' in line_out:
            label_out[iline] = 'SPACECRAFT_EVENT_STOP_TIME           = '+ spm2date(year,doy,scet[-1]) +' ' 
        iline += 1

    outlblfile = outfile[0:-4]+'.LBL'
    format_str = ('%s')
    with open(outlblfile,'wb') as file:
        for line in label_out:
            s = line[0:-1] + '\r\n'
            encoded_string = s.encode('utf-8')
            file.write(encoded_string)
    if not silent:
        print('wrote',outlblfile)

# create versions with psitype added
    outfile_psitype = outfile.replace('.TAB','_'+psitype+'.TAB')
    cp_cmd = 'cp '+outfile+'  '+outfile_psitype
#    print(cp_cmd)
    os.system(cp_cmd)

    outlblfile_psitype = outlblfile.replace('.LBL','_'+psitype+'.LBL')
    cp_cmd = 'cp '+outlblfile+'  '+outlblfile_psitype
#    print(cp_cmd)
    os.system(cp_cmd)

    if not silent:
        print('Wrote',os.path.basename(outfile_psitype))
        print('Wrote',os.path.basename(outlblfile_psitype))

    outdir = os.path.dirname(outfile)
    if tmp2output:
        dirname = os.path.dirname(outfile)+'/'
        outdir = dirname.replace('/tmp','/output')
        if tmp2globaloutput:
            outdir = outdir.replace('_local','')
        os.makedirs(outdir,exist_ok=True) # added 2025 Dec 24
        if not silent:
            print("Copying files to",outdir)
        infile = outfile
        inlblfile = outlblfile
        # add inversion range if requested
        infile_psitype = outfile_psitype
        inlblfile_psitype = outlblfile_psitype
        if add_inversion_range:
            outfile = outfile.replace('M_','M_'+f'{int(inversion_range[0]):06d}-{int(inversion_range[1]):06d}_')
            outlblfile = outlblfile.replace('M_','M_'+f'{int(inversion_range[0]):06d}-{int(inversion_range[1]):06d}_')
#        print('\nBefore updating version:',outfile)
        outfile_updated = update_file_version(outdir+os.path.basename(outfile))
#        print('After updating version:',os.path.basename(outfile_updated),'\n')
        cp_cmd = 'cp '+infile+' '+outdir + os.path.basename(outfile_updated)
        outfile = outfile_updated
#        print(cp_cmd)
        os.system(cp_cmd)
        if not silent:
            print(cp_cmd)

        outfile_psitype_updated = outfile_updated.replace('.TAB','_'+psitype+'.TAB')
        cp_cmd = 'cp '+infile_psitype+' '+outdir + os.path.basename(outfile_psitype_updated)
        outfile_psitype = outfile_psitype_updated
#        print('\n outfile_psitype_updated',os.path.basename(outfile_psitype_updated))
#        print(cp_cmd)
        os.system(cp_cmd)
        if not silent:
            print(cp_cmd)

#        print('\nBefore updating version',outlblfile)
        outlblfile_updated = update_file_version(outdir+os.path.basename(outlblfile))
#        print('After updating version:',os.path.basename(outlblfile_updated),'\n')
        cp_cmd = 'cp '+inlblfile+' '+outdir + os.path.basename(outlblfile_updated)
        outlblfile = outlblfile_updated
#        print(cp_cmd)
        os.system(cp_cmd)
        if not silent:
            print(cp_cmd)
           
        outlblfile_psitype_updated = outfile_psitype_updated.replace('.TAB','.LBL')
        cp_cmd = 'cp '+inlblfile_psitype+' '+outdir + os.path.basename(outlblfile_psitype_updated)
        outlblfile_psitype = outlblfile_psitype_updated
#        print(cp_cmd)
        os.system(cp_cmd)
        if not silent:
            print(cp_cmd)
          
        if not silent:
            print('\n')
    if clean:
        allTAUfiles = glob.glob(indir+'R*TAU*')
        path = indir+'old/'
        os.makedirs(path, exist_ok=True)
        for o in allTAUfiles:
            mv_cmd = 'mv '+ o+' '+path
            os.system(mv_cmd)
        print('merge_sorted_tau_files: temporary TAU files moved to '+path)
       
    return outfile,outlblfile,outfile_psitype,outlblfile_psitype,outdir
    
def min_required_dlp_res_km(res_km,res_factor):
    '''
    Return minimum required DLP resolution for a requested res_km and res_factor
    '''
    return res_km * res_factor
  
def min_required_dr_km_desired(res_km,resolution_factor):
    return res_km * resolution_factor/2

def min_valid_res_km(dlp_file,res_factor):
    '''
    Return minimum achievable resolution dr_km for a given DLP file
    '''
    # print(dlp_file)
    # print(dlp_file[dlp_file.index('DLP_')+4:dlp_file.index('DLP_')+8])
    dlp_res_km = float(dlp_file[dlp_file.index('DLP_')+4:dlp_file.index('DLP_')+8])/1000.
    # print('dlp_res_km',dlp_res_km)
    return dlp_res_km / res_factor
    
def MTR86_Fig10(figsize=(6,7),figfile = global_path_to_local_figs+'MTR86_Fig10.jpg'):
    bvals = np.logspace(-1,1,100)
    fig,ax=plt.subplots(figsize=figsize)
    ax.set_xlabel('PARAMETER b')
    ax.set_ylabel('RESOLUTION RATIO')
    ax.set_title('MTR86 Fig. 10')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(.1,10)
    ax.set_ylim(.1,25)
    ax.grid('True')

    DeltaR_W = 1.
    F = 1.
    W = 100
    Weff = W/1.65
    rhodot = 10.
    
    DeltaR_W = 2*F**2/Weff
    ratio1 = DeltaR_phi_MTR32(bvals,DeltaR_W=DeltaR_W)/DeltaR_W
    ax.plot(bvals,ratio1,color='k')
    ax.text(3.5,4,r'$\frac{\Delta R_\phi}{\Delta R_W}$',fontsize=16)
    
    Rinfty = DeltaR_inf_MTR3(bvals,Weff,F)
    ratio2 = DeltaR_phi_MTR32(bvals,DeltaR_W=DeltaR_W)/Rinfty
    ax.plot(bvals,ratio2,color='k')
    ax.text(.6,6,r'$\frac{\Delta R_\phi}{\Delta R_\infty}$',fontsize=16)

    ratio3 = DeltaR_W/Rinfty
    ax.plot(bvals,ratio3,color='k')
    ax.text(2,.5,r'$\frac{\Delta R_W}{\Delta R_\infty}$',fontsize=16)
    plt.savefig(figfile)
    print(figfile)
    plt.show()

def MTR86_Fig11(figsize=(6,6),figfile = global_path_to_local_figs+'MTR86_Fig11.jpg'):
    DeltaR_Wvals_m = np.logspace(1.2,3,1000)
    DeltaR_Wvals = DeltaR_Wvals_m /1000 # km
    # MTR86 Table I
    F_innerC = 9.2431
    rhodotC = 51.41
    F_outerA = 15.29583
    rhodotA = 83.49

    F_epsI = 1.6
    F_epsE = 2.3
    rhodoteps=8.2

    allan_dev_1sec = 5e-12 # p. 140 MTR86

    WeffC = 2*F_innerC**2/DeltaR_Wvals
    WC = 1.65*WeffC
    WeffA = 2*F_outerA**2/DeltaR_Wvals
    WA = 1.65*WeffA
    WeffepsI = 2*F_epsI**2/DeltaR_Wvals
    WepsI = 1.65*WeffepsI
    WeffepsE = 2*F_epsE**2/DeltaR_Wvals
    WepsE = 1.65*WeffepsE
    
    f0=8415e6
    bvalsC = b_MTR33(f0=f0,allan_dev_1sec=allan_dev_1sec,rhodot=rhodotC,W=WC)
    DeltaR_phi_C = DeltaR_phi_MTR32(bvalsC,DeltaR_W=DeltaR_Wvals)

    bvalsA = b_MTR33(f0=f0,allan_dev_1sec=allan_dev_1sec,rhodot=rhodotA,W=WA)
    DeltaR_phi_A = DeltaR_phi_MTR32(bvalsA,DeltaR_W=DeltaR_Wvals)

    bvalsepsI = b_MTR33(f0=f0,allan_dev_1sec=allan_dev_1sec,rhodot=rhodoteps,W=WepsI)
    DeltaR_phi_epsI = DeltaR_phi_MTR32(bvalsepsI,DeltaR_W=DeltaR_Wvals)
    
    bvalsepsE = b_MTR33(f0=f0,allan_dev_1sec=allan_dev_1sec,rhodot=rhodoteps,W=WepsE)
    DeltaR_phi_epsE = DeltaR_phi_MTR32(bvalsepsE,DeltaR_W=DeltaR_Wvals)
    
    fig,ax=plt.subplots(figsize=figsize)
    ax.set_xlabel(r'$\Delta R_W=2F^2/W_{eff}$ (m)')
    ax.set_ylabel(r'$\Delta R_\phi$ (m)')
    ax.set_title('MTR86 Fig. 11')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(10,1000)
    ax.set_ylim(10,1000)
    ax.grid('True')

    ax.plot(DeltaR_Wvals_m, DeltaR_Wvals_m,linestyle='dashed',linewidth=1,color='k')
    reslim=0.2
    LsolidA = np.where(DeltaR_phi_A>reslim)
    LdashedA = np.where(DeltaR_phi_A<=reslim)
    LsolidC = np.where(DeltaR_phi_C>reslim)
    LdashedC = np.where(DeltaR_phi_C<=reslim)
    LsolidepsI = np.where(DeltaR_phi_epsI>30/1000)
    LdashedepsI = np.where(DeltaR_phi_epsI<=30/1000)
    LsolidepsE = np.where(DeltaR_phi_epsE>45/1000)
    LdashedepsE = np.where(DeltaR_phi_epsE<=45/1000)
    p=ax.plot(DeltaR_Wvals_m[LsolidA], DeltaR_phi_A[LsolidA]*1000,label='Outer A ring',color='k')
    ax.plot(DeltaR_Wvals_m[LdashedA], DeltaR_phi_A[LdashedA]*1000,linestyle='dashed',color=p[0].get_color())
    p=ax.plot(DeltaR_Wvals_m[LsolidC], DeltaR_phi_C[LsolidC]*1000,label='Inner C ring',color='k')
    ax.plot(DeltaR_Wvals_m[LdashedC], DeltaR_phi_C[LdashedC]*1000,linestyle='dashed',color=p[0].get_color())

    ax.text(20,120,'OUTER A',rotation=10)
    ax.text(20, 80,'INNER C',rotation=12)
    ax.text(12,70,'SATURN',rotation=90)
    
    p=ax.plot(DeltaR_Wvals_m[LsolidepsI], DeltaR_phi_epsI[LsolidepsI]*1000,label=r'ENTER - $\epsilon$',color='k')
    ax.plot(DeltaR_Wvals_m[LdashedepsI], DeltaR_phi_epsI[LdashedepsI]*1000,linestyle='dashed',color=p[0].get_color())
    p=ax.plot(DeltaR_Wvals_m[LsolidepsE], DeltaR_phi_epsE[LsolidepsE]*1000,label=r'EXIT - $\epsilon$',color='k')
    ax.plot(DeltaR_Wvals_m[LdashedepsE], DeltaR_phi_epsE[LdashedepsE]*1000,linestyle='dashed',color=p[0].get_color())

    ax.text(17,37,r'EXIT$-\epsilon$',rotation=30)
    ax.text(17,29,r'EXIT$-\epsilon$',rotation=35)
    ax.text(12,22,'URANUS',rotation=90)
    
    #ax.legend()
    plt.savefig(figfile)
    print(figfile)
    plt.show()
    
def MTR86_Fig11_Cassini(figsize=(6,6),title='MTR86 Fig. 11 Cassini Rev137E',figfile = global_path_to_local_figs+'MTR86_Fig11_Cassini.jpg'):
    DeltaR_Wvals_m = np.logspace(1.2,3,1000)
    DeltaR_Wvals = DeltaR_Wvals_m /1000 # km
    # MTR86 Table I
    # F_innerC = 9.2431
    # rhodotC = 51.41
    # F_outerA = 15.29583
    # rhodotA = 83.49

    # # Rev007E RSS_2005_123_X43_E_GEO.TAB
    # F_innerC =2.93
    # rhodotC = 13.
    # F_outerA = 2.62
    # rhodotA = 13.07

    # Rev137E RSS_2010_245_X25_E_GEO.TAB
    F_innerC = 25.55
    rhodotC = 5.74
    F_outerA = 16.79
    rhodotA = 7.24
    
    allan_dev_1sec = 2e-13 # John Armstrong email

    WeffC = 2*F_innerC**2/DeltaR_Wvals
    WC = 1.65*WeffC
    WeffA = 2*F_outerA**2/DeltaR_Wvals
    WA = 1.65*WeffA
    
    bvalsC = b_MTR33(allan_dev_1sec=allan_dev_1sec,rhodot=rhodotC,W=WC)
    DeltaR_phi_C = DeltaR_phi_MTR32(bvalsC,DeltaR_W=DeltaR_Wvals)

    bvalsA = b_MTR33(allan_dev_1sec=allan_dev_1sec,rhodot=rhodotA,W=WA)
    DeltaR_phi_A = DeltaR_phi_MTR32(bvalsA,DeltaR_W=DeltaR_Wvals)
 
    fig,ax=plt.subplots(figsize=figsize)
    ax.set_xlabel(r'$\Delta R_W=2F^2/W_{eff}$ (m)')
    ax.set_ylabel(r'$\Delta R_\phi$ (m)')
    ax.set_title(title)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(10,1000)
    ax.set_ylim(10,1000)
    ax.grid('True')

    ax.plot(DeltaR_Wvals_m, DeltaR_Wvals_m,linestyle='dashed',linewidth=1,color='k')
    reslim=0.1
    LsolidA = np.where(DeltaR_phi_A>reslim)
    LdashedA = np.where(DeltaR_phi_A<=reslim)
    LsolidC = np.where(DeltaR_phi_C>reslim)
    LdashedC = np.where(DeltaR_phi_C<=reslim)
    p=ax.plot(DeltaR_Wvals_m[LsolidA], DeltaR_phi_A[LsolidA]*1000,label='Outer A ring')
    ax.plot(DeltaR_Wvals_m[LdashedA], DeltaR_phi_A[LdashedA]*1000,linestyle='dashed',color=p[0].get_color())
    p=ax.plot(DeltaR_Wvals_m[LsolidC], DeltaR_phi_C[LsolidC]*1000,label='Inner C ring')
    ax.plot(DeltaR_Wvals_m[LdashedC], DeltaR_phi_C[LdashedC]*1000,linestyle='dashed',color=p[0].get_color())

    # ax.text(20,120,'OUTER A',rotation=10)
    # ax.text(20, 80,'INNER C',rotation=12)
    # ax.text(12,70,'Cassini',rotation=90)
    
    ax.legend()
    plt.savefig(figfile)
    print(figfile)
    plt.show()

# particle size distribution
def n(a,q):
    return a**-q
def particle_size_models(w_K = 9e6,w_X = 3.6e7,w_S = 13.e7,q_vals = [2.8,3.0,3.2,3.4],
                         m_vals = [1.78+0j],a_min_vals = np.logspace(6,9,301),a_max_vals = [3e9,1e10]):
    '''
    ~~ Compute Model Grids ~~
    '''
    
    # storage lists for plotting
    tau_ratio_KX = [] # (tau_K - tau_X)/tau_x
    tau_ratio_XS = [] # (tau_X - tau_S)/tau_X
    labels = []  # labels for plot legend
    colors = []  # colors for plot lines
    markers = [] # markers for plot symbols
    lines = []   # line style for plot lines
    
    # compute efficiencies for difference particle compositions
    for m,ms,ls in zip(m_vals,['o','s'],['-','--']):
    
        ## Get differential Mie extinction efficiencies for full range of radii
        # set PyMieScatt keywargs
        drange = (2*a_min_vals[0],2*a_max_vals[-1])
        kwargs = {'nMedium':1.0,'diameterRange':drange,'nd':1000,'logD':True}
        # use PyMieScatt to compute Mie extinction efficiencies at three wavelengths
        d,Q_K = ps.MieQ_withDiameterRange(m,w_K,**kwargs)[0:2]
        d,Q_X = ps.MieQ_withDiameterRange(m,w_X,**kwargs)[0:2]
        d,Q_S = ps.MieQ_withDiameterRange(m,w_S,**kwargs)[0:2]
        # differential extinction efficiencies from extinction efficiencies
        dQ_KX = Q_K - Q_X
        dQ_XS = Q_X - Q_S
        # particle radii from diameters
        a_vals = d/2
    
        # compute integrals for different radii ranges
        for a_max,c in zip(a_max_vals,['0.0','0.6']):
    
            # compute integrals for different power law slopes
            for qi in q_vals:
    
                # new empty lists for storage
                Tau_X = []
                DeltaTau_KX = []
                DeltaTau_XS = []
    
                # compute dtau/tau for a range of minimum particle radii
                for a_min in a_min_vals:
    
                    # clipping mask array based on particle radius
                    rclip = np.where((a_vals>=a_min)&(a_vals<=a_max))[0]
                    # compute tau integral with Simpson's rule and store
                    Tau_X += [simpson(tau_int(a_vals[rclip],Q_X[rclip],q=qi),x=a_vals[rclip])]
                    # compute Delta tau integrals with Simpson's rule and store in lists
                    DeltaTau_KX += [simpson(dtau_int(a_vals[rclip],dQ_KX[rclip],q=qi),x=a_vals[rclip])]
                    DeltaTau_XS += [simpson(dtau_int(a_vals[rclip],dQ_XS[rclip],q=qi),x=a_vals[rclip])]
    
                # convert lists to numpy arrays
                Tau_X = np.array(Tau_X)
                DeltaTau_KX = np.array(DeltaTau_KX)
                DeltaTau_XS = np.array(DeltaTau_XS)
    
                # store arrays for plotting
                tau_ratio_KX += [DeltaTau_KX/Tau_X]
                tau_ratio_XS += [DeltaTau_XS/Tau_X]
    
                # store plotting parameters for later
                markers += [ms]
                lines += [ls]
                colors += [c]
                labels += [r'$a_{max}=$'+str(a_max/1e9)+' m']#[r'$\overline{m}=$'+str(np.real(m))+r', $a_{max}=$'+str(a_max/1e9)+' m']
    return q_vals,m_vals,tau_ratio_KX,tau_ratio_XS,markers,lines,colors,labels
    
def pick_SC_kernel(CORSS_TABfile,kernels_list = global_path_to_tables+'e2e_kernels.ker',verbose=False):
    '''
    Pick a Cassini S/C ephemeris file from tabulated list, parsing CORSS_TABfile for date of event
    '''
    kernel_found = "NO KERNEL FOUND" # default response if no kernel found
    index_start = CORSS_TABfile.index('_RSS_')
    frag = CORSS_TABfile[index_start+5:]
    year_start = frag[2:4]
    doy_start = frag[5:8]
    yeardoy_start = int(year_start + doy_start)

    if verbose:
        print(CORSS_TABfile)
    
    with open(kernels_list,'r') as f:
        kernels_all  = f.read().splitlines() # strips \n at end of each line.
    for kernel in kernels_all:
        try: # skip over lines in kernels_all that aren't Cassini S/C ephemerides
            index= kernel.index('SCPSE')
            # get date range covered by this kernel
            yrddds= int(kernel[index+6:index+11])
            yrddde= int(kernel[index+12:index+17])
            if (yeardoy_start >= yrddds) and (yeardoy_start <= yrddde):
                #print('debug pick_SC_kernel')
                #print(kernel)
                #print('yrddds,yrddde',yrddds,yrddde)
                index = kernel.index('$A/') +len('$A/')
                kernel_found = kernel[index:index+len('050623R_SCPSE_05132_05150.bsp')]
                break
        except: # this line does not contain a Cassini S/C ephemeris
            continue
    return 'naif/CASSINI/kernels/spk/' + kernel_found
    
def pick_RSR_file(_16kHz=False,Rev='Rev007',direction='E',DSN='X43',
                  local_path_to_tables=global_path_to_tables,force=False,download=False,silent=True):
    '''
    Grab an RSR file for a given Rev, direction, DSN (band/dsn number)
    '''
    #print("TP1 - pick_RSR_file: _16kHz=",_16kHz)
    if _16kHz:
        file_list = local_path_to_tables+'rsr_16kHz_files_before_USO_failure.txt'
    else:
        file_list = local_path_to_tables+'rsr_all_files_before_USO_failure.txt'
    with open(file_list,'r') as f:
        RSR_files = f.read().splitlines()
    if not _16kHz:
        rsr_files =[]
        for rsr_file in RSR_files:
            if rsr_file[-1]=='1': # suffix 1 means 1 kHz
                rsr_files.append(rsr_file)
        RSR_files = rsr_files

    CORSS_TABfile = get_CORSS_8001_TABfiles(CORSS_8001_all_filepaths=None,
            Rev=Rev,DSN=DSN,direction=direction,silent=True)
    if not silent:
        print("CORSS_TABfile",CORSS_TABfile)
    
    rsr_file = 'No matching RSR file found for this event'
    if len(CORSS_TABfile) != 0:
        event= os.path.basename(CORSS_TABfile)
        string = direction + event[4:8]+event[9:12] # string to be found in desired RSR
        for rsr_file in RSR_files:
            if (string.upper() in rsr_file.upper()) and (DSN.upper() in rsr_file.upper()):
                if not silent:
                    print('Using RSR file',rsr_file)
                get_RSRfiles_from_web(rsr_file,force=force,silent=True) # make sure we have it
                return rsr_file
# this can fail for chord occ if searching for egress, since RSR file for chord occ is ingress
    if direction =='E' and CORSS_TABfile.index('CE') >0: # 'CE' indicates chord egress      
        string = 'I' + event[4:8]+event[9:12] # search for ingress instead
        for rsr_file in RSR_files:
            if (string.upper() in rsr_file.upper()) and (DSN.upper() in rsr_file.upper()):
                if not silent:
                    print('Using RSR file',rsr_file)
                get_RSRfiles_from_web(rsr_file,silent=True) # make sure we have it
                return rsr_file
# it can also fail if 1 kHz file requested but not found
    if (len(CORSS_TABfile) != 0) and not _16kHz:
        event= os.path.basename(CORSS_TABfile) #copied from above for clarity
        string = direction + event[4:8]+event[9:12] # string to be found in desired RSR
        file_list = local_path_to_tables+'rsr_16kHz_files_before_USO_failure.txt'
        with open(file_list,'r') as f:
            RSR_files = f.read().splitlines()
        rsr_files =[]
        for rsr_file in RSR_files:
            if (string.upper() in rsr_file.upper()) and (DSN.upper() in rsr_file.upper()):
                print('No 1 kHz RSR file found, so using 16 kHz file',rsr_file)
                get_RSRfiles_from_web(rsr_file,force=force,silent=True) # make sure we have it
                return rsr_file
    if string == 'E2008102' and DSN=='K47':
        rsr_file = 'co-s-rss-1-sroc7-v10/cors_0225/sroc7_102/rsr/s39sroi2008102_0705nnnk47rd.2a2'
#        print('DEBUG: setting rsr_file to ',rsr_file)
        return rsr_file
    if string == 'E2008232' and DSN=='K47':
        rsr_file = 'co-s-rss-1-sroc8-v10/cors_0249/sroc8_232/rsr/s43sroi2008232_0517nnnk47rd.2a2'
#        print('DEBUG: setting rsr_file to ',rsr_file)
        return rsr_file
    print('RSR file not found for',Rev,direction,DSN)
    return None

    
def plot_comparisons(program,results,show=False,y_offset=0.1,dy_offset=0.,plot_phase=False,tau=None,
                     adjust_phase=False,xlim=(None,None),ylim_tau=(None,None),ylim_phase=(None,None),
                     plot_tau_threshold=True,
                     CORSS_8001_filepath=None,include_CORSS_8001=False,figsize=(8,5),fontsize_legend = 'x-small',
                     ncol_legend=2,loc_legend='lower left',
                     local_path_to_data=global_path_to_local_data,path_to_figs=global_path_to_local_figs,
                     title=None,save_figfile=True,this_figfilename=None,verbose=False):
    """
    loop over results dictionaries and plot tau(r), including CORSS_8001 comparison if supplied
    results dict contains required information about DLP, TAB and CORSS_8001 and reconstruction info
    
    program: string - name of calling program, used as prefix for figure filename
    results: dict of results returned from 
    y_offset applies this offset to plotted variable 
    dy_offset applies this offset to succesive plots in results array 
    plot_phase plots phase instead of optical depth
    adjust_phase aligns phase of first point of plotted reconstruction with corresponding phase from CORSS_8001,
    since an arbitrary constant phase does not affect the reconstruction and this permits an easier comparison
    of the post-reconstruction phases in CORSS_8001 and rss_ringoccs results

    path_to_figs: is created if not already present

    figfilename has unique data string unless thisfigfilename keyword is set
    """
    os.makedirs(path_to_figs, exist_ok=True) # check to see if required figs directory exists
        
    if type(results) == type({'a':0}):
        results= [results] # ensure it is a list of results, needed for loop below even if only one results supplied
    try:
        rr = results[0] 
    except:
        if verbose:
            print("Error in results in plot_comparisons") # happens when requested data range not present
        return 'NOFIGFILE'
    dlp_file = rr['dlp_file']
    psitype = rr['psitype']
    inversion_range = rr['inversion_range']

    # string = os.path.basename(dlp_file)
    # index = string.index('_DLP')
    # event_name = string[0:index]
    Rev = rr['Rev']
    direction = rr['direction']
    DSN = rr['DSN']
    plot_title = 'RSS_'+Rev[3:]+'_'+DSN+direction
    if type(title) == type('a string'):
        plot_title = plot_title +'   '+title

    if plot_phase:
        index = 7 
        ylabel= 'Phase (deg)'
        if adjust_phase:
            ylabel = 'Adjusted '+ ylabel
        fig_suffix = '_phase_'
    else:
        index = 6
        ylabel = 'Normal optical depth'
        fig_suffix = '_tau_'
    plt.figure(figsize=figsize)
        
    if title != None:
        plt.title(title)
    plt.grid(axis='x',linestyle='--')
    plt.xlabel('Radius (km)')
    plt.ylabel('Optical Depth')
    lw_CORSS = 2
    lw_RSS_RINGOCCS=1
    for i,rr in enumerate(results):
        Rev = rr['Rev']
        direction = rr['direction']
        DSN = rr['DSN']
        psitype = rr['psitype']
        res_km = rr['res_km']
        res_factor = rr['res_factor']
        wtype = rr['wtype']
        band = rr['band']
        label='RSS'+Rev[3:]+'_'+band+DSN+direction+' '+psitype+' '+str(res_km)+ ' km ('+str(res_factor)+') '+wtype
        if rr['tau_file'] == None:
            rkm = rr['rkm']
            tau = rr['tau']
            tau_threshold = rr['tau_threshold']
        else:
            tau_file = rr['tau_file']
            d = np.genfromtxt(tau_file,delimiter=',',skip_header=1)
            rkm = d[:,0] # get radius from tau fie
            tau = d[:,index] + y_offset
            tau_threshold = d[:,8]

       # print('*** in plot_comparisons: len(rkm),len(tau_threshold):',len(rkm),len(tau_threshold))
        
        yvals = tau + y_offset
        if i == 0 and rr['CORSS_8001_filepath'] != None and include_CORSS_8001 == True:
            local_CORSS_8001_TAUfile = get_CORSS_8001_file(rr['CORSS_8001_filepath'],\
                 local_path_to_data = local_path_to_data,silent=True)
#            print('DEBUG: local_CORSS_8001_TAUfile=',local_CORSS_8001_TAUfile)
            d = np.genfromtxt(local_CORSS_8001_TAUfile,delimiter=',')
            rPDS3 = d[:,0] # radius
            yPDS3 = d[:,index] # optical depth or phase
            plt.plot(rPDS3,yPDS3,linewidth=lw_CORSS,label='CORSS_8001 '+'RSS_'+Rev[3:]+'_'+DSN+direction)
        if plot_phase:
            if adjust_phase:
                LPDS3 = min(np.where(rPDS3 >= min(inversion_range)))[0]
                L = min(np.where(rkm >= min(inversion_range)))[0]
                dymatch = yPDS3[LPDS3]-yvals[L]
            else:
                dymatch = 0
            yvals = (yvals+dymatch + 360+180)%360 - 180
        if i==0:
            p=plt.plot(rkm,yvals,label=label,linewidth=lw_RSS_RINGOCCS)               
        else:
            p=plt.plot(rkm,yvals+i*dy_offset,linewidth=0.75,label=label)
        if not plot_phase and plot_tau_threshold:
            plt.plot(rkm,tau_threshold+i*dy_offset,linestyle='dashed',color=p[0].get_color())
    if plot_phase and not adjust_phase:
        ymin_plot=-180
        ymax_plot = 180
    else:
        ymax = np.max(yvals+i*dy_offset)
        ymin = np.min(yvals+i*dy_offset)
        if ymax >0:
            ymax_plot = 1.2*ymax
        else:
            ymax_plot = 0.8*ymax
        if ymin >0:
            ymin_plot = 0.8*ymin
        else:
            ymin_plot = 1.2*ymin
        if not plot_phase:
            ymin_plot=-.05 # for optical depth
    ylim = (ymin_plot,ymax_plot)
    if ylabel == 'Normal optical depth' and ylim_tau != (None,None):
        ylim = ylim_tau
    if fig_suffix == '_phase_' and ylim_phase != (None,None):
        ylim = ylim_phase
    plt.ylim(ylim)
    plt.legend(fontsize=fontsize_legend,loc=loc_legend,ncol=ncol_legend)
    if xlim == (None,None):
        plt.xlim(inversion_range)
    else:
        plt.xlim(xlim)
    plt.xlabel('Radius (km)')
    plt.ylabel(ylabel)

    # timestamped figure filename 
    if this_figfilename == None:
            figfile = path_to_figs+program+fig_suffix+time.strftime("%Y-%m-%d_%H:%M:%S")+'.png' 
    else:
        figfile = path_to_figs + this_figfilename
    if save_figfile:
        plt.savefig(figfile)
        print('Figure saved as',figfile)
    if show:
        plt.show()
    return figfile
    
def printnow(CRbefore=True,CRafter=True):
    '''
    Print current time in nice format.
    '''
    current_time = datetime.datetime.now()
    if CRbefore:
        prefix = '\n'
    else:
        prefix = ''
    if CRafter:
        suffix = '\n'
    else:
        suffix = ''
    print(prefix+"Current time:", current_time.strftime("%H:%M:%S")+suffix)

def print_featurelist(d):
    '''
    Print a nicely formatted list of features
    '''
    print('Features:')
    print('index   name            a_res   dr_min  dr_max')
    for index,name,a_res,dr_min,dr_max in zip(d['index'],d['name'],d['a_res'],d['dr_min'],d['dr_max']):
        print(f'{index:3d}  {name:15}  {a_res:9.1f} {dr_min:4.0f}    {dr_max:4.0f}')                                 

def radius_correction_pole(tau_file,fit_number=1,verbose=False,NMAX=500):
    '''
    Compute radius correction due to difference in pole between nominal kernel and NCFIV fits 1 or 7
    '''
    rev_info = get_rev_info_from_tau(tau_file)
    UTCdate = rev_info['year']+'-'+rev_info['doy']+'T00:00'
    ETdate = spice.str2et(UTCdate)
    rho_km_vals_taufile,spm_OET_taufile = np.loadtxt(tau_file,delimiter=',',usecols=[0,9]).T
    et_vals = ETdate + spm_OET_taufile
    planet ='Saturn'
    spacecraft = 'Cassini'
    dsn = 'DSS-'+rev_info['dsn'] 
    Npts = len(rho_km_vals_taufile)
    if Npts > NMAX:
        et_vals_ = np.linspace(et_vals[0],et_vals[-1],NMAX)
    else:
        et_vals_ = et_vals
        
# recalculate rho_km_vals - these should be very close to rkm_taufile values
    rho_km_vals_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)

    if Npts > NMAX:
        frho_of_et = interpolate.interp1d(et_vals,rho_km_vals_taufile,
                                          bounds_error=False,fill_value='extrapolate')
        rho_km_vals_calc = frho_of_et(et_vals_)
    else:
        rho_km_vals_calc = rho_km_vals_
    if verbose:
        plt.plot(rho_km_vals_calc,rho_km_vals_calc-rho_km_vals_taufile,label='recalc - PDS')
        plt.legend()
        plt.show()

# save input pole direction
    BODY699_POLE_RA = spice.gdpool('BODY699_POLE_RA',0,3)
    BODY699_POLE_DEC = spice.gdpool('BODY699_POLE_DEC',0,3)

# NCFIV pole fits
    assert (fit_number==1 or fit_number==7),'radius_correction_pole: illegal fit_number '+str(fit_number)
    if fit_number == 1:
        BODY699_POLE_RA_fit = np.array([40.579414,-.03497,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537218,-.00324,0]) # at epoch 2008 Jan 1 12:00 UTC
    elif fit_number == 7:
        BODY699_POLE_RA_fit = np.array([40.579425,-.03062,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537202,-.00461,0]) # at epoch 2008 Jan 1 12:00 UTC
# correct these for J2000
    epoch = '2008 Jan 1 12:00'
    det_epoch_cy = spice.str2et(epoch)/36525/86400
# replace pole with fitted pole from NCFIV
    BODY699_POLE_RA_fit[0] -= BODY699_POLE_RA_fit[1]*det_epoch_cy
    BODY699_POLE_DEC_fit[0] -= BODY699_POLE_DEC_fit[1]*det_epoch_cy
# update the pole direction in the kernel pool
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA_fit)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC_fit)
# calculate revised radius scale
    rho_km_vals_corr_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)
    drho_km_vals_corr_ = rho_km_vals_corr_ - rho_km_vals_calc
    if Npts > NMAX:
        fdrho_et = interpolate.interp1d(et_vals_,drho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        drho_km_vals_corr = fdrho_et(et_vals)
        frho_corr_et = interpolate.interp1d(et_vals_,rho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        rho_km_vals_corr =frho_corr_et(et_vals)
    else:
        drho_km_vals_corr = drho_km_vals_corr_
        rho_km_vals_corr = rho_km_vals_corr_
        
# restore original pole direction
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC)
    
    return drho_km_vals_corr,rho_km_vals_calc,rho_km_vals_corr
    
def radius_correction_pole_from_tau_inst(dlp_file,tau_inst,fit_number=1,verbose=False,NMAX=500):
    '''
    Compute radius correction due to difference in pole between nominal kernel and NCFIV fits 1 or 7
    '''
    rev_info = get_rev_info_from_dlp(dlp_file)
    UTCdate = rev_info['year']+'-'+rev_info['doy']+'T00:00'
    ETdate = spice.str2et(UTCdate)
    rho_km_vals_taufile = tau_inst.rho_km_vals
    spm_OET_taufile = tau_inst.t_oet_spm_vals
    et_vals = ETdate + spm_OET_taufile
    planet ='Saturn'
    spacecraft = 'Cassini'
    dsn = 'DSS-'+rev_info['dsn'] 
    Npts = len(rho_km_vals_taufile)
    if Npts > NMAX:
        et_vals_ = np.linspace(et_vals[0],et_vals[-1],NMAX)
    else:
        et_vals_ = et_vals
        
# recalculate rho_km_vals - these should be very close to rkm_taufile values
    rho_km_vals_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)

    if Npts > NMAX:
        frho_of_et = interpolate.interp1d(et_vals,rho_km_vals_taufile,
                                          bounds_error=False,fill_value='extrapolate')
        rho_km_vals_calc = frho_of_et(et_vals_)
    else:
        rho_km_vals_calc = rho_km_vals_
    if verbose:
        plt.plot(rho_km_vals_calc,rho_km_vals_calc-rho_km_vals_taufile,label='recalc - PDS')
        plt.legend()
        plt.show()

# save input pole direction
    BODY699_POLE_RA = spice.gdpool('BODY699_POLE_RA',0,3)
    BODY699_POLE_DEC = spice.gdpool('BODY699_POLE_DEC',0,3)

# NCFIV pole fits
    assert (fit_number==1 or fit_number==7),'radius_correction_pole: illegal fit_number '+str(fit_number)
    if fit_number == 1:
        BODY699_POLE_RA_fit = np.array([40.579414,-.03497,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537218,-.00324,0]) # at epoch 2008 Jan 1 12:00 UTC
    elif fit_number == 7:
        BODY699_POLE_RA_fit = np.array([40.579425,-.03062,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537202,-.00461,0]) # at epoch 2008 Jan 1 12:00 UTC
# correct these for J2000
    epoch = '2008 Jan 1 12:00'
    det_epoch_cy = spice.str2et(epoch)/36525/86400
# replace pole with fitted pole from NCFIV
    BODY699_POLE_RA_fit[0] -= BODY699_POLE_RA_fit[1]*det_epoch_cy
    BODY699_POLE_DEC_fit[0] -= BODY699_POLE_DEC_fit[1]*det_epoch_cy
# update the pole direction in the kernel pool
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA_fit)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC_fit)
# calculate revised radius scale
    rho_km_vals_corr_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)
    drho_km_vals_corr_ = rho_km_vals_corr_ - rho_km_vals_calc
    if Npts > NMAX:
        fdrho_et = interpolate.interp1d(et_vals_,drho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        drho_km_vals_corr = fdrho_et(et_vals)
        frho_corr_et = interpolate.interp1d(et_vals_,rho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        rho_km_vals_corr =frho_corr_et(et_vals)
    else:
        drho_km_vals_corr = drho_km_vals_corr_
        rho_km_vals_corr = rho_km_vals_corr_
        
# restore original pole direction
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC)
    
    return drho_km_vals_corr,rho_km_vals_calc,rho_km_vals_corr
    
def radius_correction_trajectory(tau_file,fit_number=1,NMAX=1000,verbose=False):
    '''
    Compute radius correction due to trajectory error from NCFIV fits 1 or 7.
    If more than NMAX points in tau_file, use interpolation to speed up calculation.
    '''
    rev_info = get_rev_info_from_tau(tau_file)
    UTCdate = rev_info['year']+'-'+rev_info['doy']+'T00:00'
    ETdate = spice.str2et(UTCdate)
    rho_km_vals_taufile,spm_OET_taufile = np.loadtxt(tau_file,delimiter=',',usecols=[0,9]).T
    et_vals = ETdate + spm_OET_taufile
    planet ='Saturn'
    spacecraft = 'Cassini'
    dsn = 'DSS-'+rev_info['dsn'] 
    
    Npts = len(rho_km_vals_taufile)
    if Npts > NMAX:
        et_vals_ = np.linspace(et_vals[0],et_vals[-1],NMAX)
    else:
        et_vals_ = et_vals
        
# recalculate rho_km_vals - these should be very close to rkm_taufile values
    rho_km_vals_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)

    if Npts > NMAX:
        frho_of_et = interpolate.interp1d(et_vals_,rho_km_vals_,
                                          bounds_error=False,fill_value='extrapolate')
        rho_km_vals_calc = frho_of_et(et_vals)
    else:
        rho_km_vals_calc = rho_km_vals_
    if verbose:
        plt.plot(rho_km_vals_calc,rho_km_vals_calc-rho_km_vals_taufile,label='recalc - PDS')
        plt.legend()
        plt.show()

# save input pole direction
    BODY699_POLE_RA = spice.gdpool('BODY699_POLE_RA',0,3)
    BODY699_POLE_DEC = spice.gdpool('BODY699_POLE_DEC',0,3)

# NCFIV pole fits
    assert (fit_number==1 or fit_number==7),'radius_correction_pole: illegal fit_number '+str(fit_number)
    if fit_number == 1:
        BODY699_POLE_RA_fit = np.array([40.579414,-.03497,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537218,-.00324,0]) # at epoch 2008 Jan 1 12:00 UTC
    elif fit_number == 7:
        BODY699_POLE_RA_fit = np.array([40.579425,-.03062,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537202,-.00461,0]) # at epoch 2008 Jan 1 12:00 UTC
# correct these for J2000
    epoch = '2008 Jan 1 12:00'
    det_epoch_cy = spice.str2et(epoch)/36525/86400
# replace pole with fitted pole from NCFIV
    BODY699_POLE_RA_fit[0] -= BODY699_POLE_RA_fit[1]*det_epoch_cy
    BODY699_POLE_DEC_fit[0] -= BODY699_POLE_DEC_fit[1]*det_epoch_cy
# update the pole direction in the kernel pool
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA_fit)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC_fit)
# calculate revised radius scale
    rho_km_vals_corr_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)
# obtain trajectory correction coefficients from fit output file
    dt_corr,alpha,r0 = get_trajectory_correction_coefficients(rev_info,fit_number,verbose=verbose)
#    dt_corr = -0.021 # all numbers are rounded in NCFIV!
# compute radius correction
    rdot_ = np.gradient(rho_km_vals_corr_,et_vals_)
    drho_km_vals_corr_ = rdot_ * dt_corr - alpha*(rho_km_vals_corr_-r0)/1000.
    if Npts > NMAX:
        fdrho_et = interpolate.interp1d(et_vals_,drho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        drho_km_vals_corr = fdrho_et(et_vals)
        frho_corr_et = interpolate.interp1d(et_vals_,rho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        rho_km_vals_corr =frho_corr_et(et_vals)
    else:
        drho_km_vals_corr = drho_km_vals_corr_
        rho_km_vals_corr = rho_km_vals_corr_
        
# restore original pole direction
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC)
    
    return drho_km_vals_corr,rho_km_vals_calc,rho_km_vals_corr
    
def radius_correction_trajectory_from_tau_inst(dlp_file,tau_inst,fit_number=1,NMAX=1000,verbose=False):
    '''
    Compute radius correction due to trajectory error from NCFIV fits 1 or 7
    '''

    rev_info = get_rev_info_from_dlp(dlp_file)

    UTCdate = rev_info['year']+'-'+rev_info['doy']+'T00:00'
    ETdate = spice.str2et(UTCdate)
    rho_km_vals_taufile = tau_inst.rho_km_vals
    spm_OET_taufile = tau_inst.t_oet_spm_vals
    et_vals = ETdate + spm_OET_taufile
    planet ='Saturn'
    spacecraft = 'Cassini'
    dsn = 'DSS-'+rev_info['dsn'] 
# the following code copied from   radius_correction_trajectory()  
    Npts = len(rho_km_vals_taufile)
    if Npts > NMAX:
        et_vals_ = np.linspace(et_vals[0],et_vals[-1],NMAX)
    else:
        et_vals_ = et_vals
        
# recalculate rho_km_vals - these should be very close to rkm_taufile values
    rho_km_vals_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)

    if Npts > NMAX:
        frho_of_et = interpolate.interp1d(et_vals_,rho_km_vals_,
                                          bounds_error=False,fill_value='extrapolate')
        rho_km_vals_calc = frho_of_et(et_vals)
    else:
        rho_km_vals_calc = rho_km_vals_
    if verbose:
        plt.plot(rho_km_vals_calc,rho_km_vals_calc-rho_km_vals_taufile,label='recalc - PDS')
        plt.legend()
        plt.show()

# save input pole direction
    BODY699_POLE_RA = spice.gdpool('BODY699_POLE_RA',0,3)
    BODY699_POLE_DEC = spice.gdpool('BODY699_POLE_DEC',0,3)

# NCFIV pole fits
    assert (fit_number==1 or fit_number==7),'radius_correction_pole: illegal fit_number '+str(fit_number)
    if fit_number == 1:
        BODY699_POLE_RA_fit = np.array([40.579414,-.03497,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537218,-.00324,0]) # at epoch 2008 Jan 1 12:00 UTC
    elif fit_number == 7:
        BODY699_POLE_RA_fit = np.array([40.579425,-.03062,0]) # at epoch 2008 Jan 1 12:00 UTC
        BODY699_POLE_DEC_fit = np.array([83.537202,-.00461,0]) # at epoch 2008 Jan 1 12:00 UTC
# correct these for J2000
    epoch = '2008 Jan 1 12:00'
    det_epoch_cy = spice.str2et(epoch)/36525/86400
# replace pole with fitted pole from NCFIV
    BODY699_POLE_RA_fit[0] -= BODY699_POLE_RA_fit[1]*det_epoch_cy
    BODY699_POLE_DEC_fit[0] -= BODY699_POLE_DEC_fit[1]*det_epoch_cy
# update the pole direction in the kernel pool
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA_fit)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC_fit)
# calculate revised radius scale
    rho_km_vals_corr_ = calc_rho_km(et_vals_, planet, spacecraft, dsn)
# obtain trajectory correction coefficients from fit output file
    dt_corr,alpha,r0 = get_trajectory_correction_coefficients(rev_info,fit_number,verbose=verbose)
#    dt_corr = -0.021 # all numbers are rounded in NCFIV!
# compute radius correction
    rdot_ = np.gradient(rho_km_vals_corr_,et_vals_)
    drho_km_vals_corr_ = rdot_ * dt_corr - alpha*(rho_km_vals_corr_-r0)/1000.
    if Npts > NMAX:
        fdrho_et = interpolate.interp1d(et_vals_,drho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        drho_km_vals_corr = fdrho_et(et_vals)
        frho_corr_et = interpolate.interp1d(et_vals_,rho_km_vals_corr_,
                                              bounds_error=False,fill_value='extrapolate')
        rho_km_vals_corr =frho_corr_et(et_vals)
    else:
        drho_km_vals_corr = drho_km_vals_corr_
        rho_km_vals_corr = rho_km_vals_corr_
        
# restore original pole direction
    spice.pdpool('BODY699_POLE_RA',BODY699_POLE_RA)
    spice.pdpool('BODY699_POLE_DEC',BODY699_POLE_DEC)
    
    return drho_km_vals_corr,rho_km_vals_calc,rho_km_vals_corr
    

def rkm_tau_from_TAUfile(path_to_tau_file):
    data = np.loadtxt(path_to_tau_file,delimiter=',')
    rkm = data[:,0]
    tau = data[:,6]
    return rkm, tau

def rkm_tau_drs_tanBeff_from_TAUfile(path_to_tau_file):
    data = np.loadtxt(path_to_tau_file,delimiter=',')
    rkm = data[:,0]
    drpole = data[:,1]
    drtraj = data[:,2]
    phi = data[:,4]
    tau = data[:,6]
    B = data[:,-1]
    tanBeff = np.tan(np.radians(B))/np.cos(np.radians(phi))
    return rkm, tau, drpole,drtraj,tanBeff
    
def rkm_tau_tanBeff_from_TAUfile(path_to_tau_file):
    data = np.loadtxt(path_to_tau_file,delimiter=',')
    rkm = data[:,0]
    phi = data[:,4]
    tau = data[:,6]
    B = data[:,-1]
    tanBeff = np.tan(np.radians(B))/np.cos(np.radians(phi))
    return rkm, tau, tanBeff
    
def RSRfiles_for_demo(_16kHz=True): 
    '''
    Grab Rev007E K34,X34,X43 16 kHz files to enable high-resolution reconstructions
    '''
    if _16kHz:
        RSRfiles = [
            'co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnk34rd.1b2',
            'co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnx34rd.1a2',
            'co-s-rss-1-sroc1-v10/cors_0105/sroc1_123/rsr/s10sroe2005123_0740nnnx43rd.2a2'
        ]
    else:
        RSRfiles = [
            'co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNK34RD.1B1',
            'co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNX34RD.1A1',
            'co-s-rss-1-sroc1-v10/cors_0727/SROC1_123/RSR/S10SROE2005123_0740NNNX43RD.2A1'
        ]
    return RSRfiles
    
def runloop(figfile=None,PRINT_TIME=False,include_history=False,verbose=False,xlim_run=(None,None),
            ylim_tau=(None,None),y_offset=0.1,dy_offset=0.,NoOp=False,
           include_CORSS_8001 = True,save_figfile = True): 
    '''
    Perform end-to-end loop over requested features. Makes use of available variables.
    '''
    tstart_all = time.time()
    results_complete = []  
    outfiles_complete = []
    if verbose:
        print('runloop: Lfeatures',Lfeatures)
    for i in Lfeatures:
        tstart_feature = time.time()
        a_res = dict_features['a_res'][i]
        dr_min = dict_features['dr_min'][i]
        dr_max = dict_features['dr_max'][i]
        name = dict_features['name'][i]
        if not NoOp and verbose:
            print('Processing',name,'...')
        rmin = a_res + dr_min
        rmax = a_res + dr_max
        inversion_range = [rmin,rmax]
        for run in runs:
            tstart_run = time.time()
    #       plt.figure(figsize=figsize)
            rev = run['rev']
            direc = run['direc']
            dlp_res = run['dlp_res']
            psitypes = run['psitypes']
            res_kms = np.array(run['res_km'])
            drho_km = run['drho_km']
            if 'xlim_run' in run:
                xlim_run = run['xlim_run']
            if 'write_tau_files' in run:
                write_tau_files = run['write_tau_files']
            else:
                write_tau_files = False
            if 'subdivide' in run and 'drange' in run:
                subdivide = run['subdivide']
                drange = run['drange']
            else:
                subdivide = False
            if 'tmp2output' in run:
                tmp2output = run['tmp2output']
            else:
                tmp2output = True
            if 'silent' in run:
                silent = run['silent']
            else:
                silent = False
            if 'clean' in run:
                clean = run['clean']
            else:
                clean = False

            search_dir = global_path_to_output+'Rev'+rev+'/Rev'+rev+'*'+direc+'/*'
            dirs = glob.glob(search_dir)
            # print('global_path_to_output:',global_path_to_output)
            # print('search_dir:',search_dir)
            # print('dirs:',dirs)
            results_all = []
           
            bands = np.array(run['bands'])
            if bands[0] == 'ALL':
                bands = bands_ALL
    
            for band in bands:
                tstart_band = time.time()
                dsns = np.array(run['dsns'])
                if dsns[0] == 'ALL':
                    dsns = dsns_ALL
                elif dsns[0] == 'OMIT 63':
                    dsns = dsns_OMIT_63
                for dsn in dsns:
                    tstart_dsn=time.time()
     #construct possible name of a directory
                    string = band + dsn+'_'+direc
#                    print('searching for directory containing '+string)
                    for this_dir in dirs:
                        if string in this_dir:
                            if verbose:
                                print('Processing',this_dir)
                            # print('now search in ',this_dir,'for GEO,CAL,DLP')
                            GEO_glob = '/*GEO*'+search_string
                            CAL_glob = '/*CAL*'+search_string
                            DLP_glob = '/*'+dlp_res+search_string
                            try:
                                geo_file = glob.glob(this_dir + GEO_glob)[0]
                            except:
                                print('Failed to get GEO file',this_dir + GEO_glob)
                                print(glob.glob(this_dir + GEO_glob))
                                geo_file = glob.glob(this_dir + GEO_glob)
                            try:
                                cal_file = glob.glob(this_dir + CAL_glob)[0]
                            except:
                                print('Failed to get CAL file',this_dir + CAL_glob)
                                print(glob.glob(this_dir + CAL_glob))
                                cal_file = glob.glob(this_dir + CAL_glob)
                            try:
                                dlp_file = glob.glob(this_dir + DLP_glob)[0]
                            except:
                                print('No matching DLP file at requested resolution',dlp_res)
                                print(this_dir + DLP_glob)
                                dlp_res = 'DLP_*???M'
                                DLP_glob_all = '/*'+dlp_res+'_'+search_string
                                dlp_files = glob.glob(this_dir + DLP_glob_all)
                                if verbose:
                                    print('found:')
                                    for dlp_file in dlp_files:
                                        print(dlp_file)
                                try:
                                    dlp_file = dlp_files[0] # use highest resolution one in list
                                except:
                                    print('No DLP file found. Likely to be an incomplete directory, skipping....')
                                    return outfiles_complete # perhaps some of these have been successfully run
                            print('\nUsing',os.path.basename(dlp_file))
                            string = os.path.basename(dlp_file)
                            index = string.index('DLP')
                            dlp_res_used = string[index:index+9]
                                    
                        # get radial range of dlp_file to see if reqested range is present
                            dlp_contents = np.loadtxt(dlp_file,delimiter=',')
                            rmin_dlp = dlp_contents[0,0]
                            rmax_dlp = dlp_contents[-1,0]
                            if (rmin_dlp > inversion_range[0]) or (rmax_dlp < inversion_range[0]):
                                print('DLP range does not span requested inversion_range')
                                print('DLP range:',rmin_dlp,rmax_dlp)
                                print('inversion_range:',inversion_range)

                            # print(os.path.basename(geo_file))
                            # print(os.path.basename(cal_file))
                            # print(os.path.basename(dlp_file))
                            # print('reading geo, cal, dlp files...')
                            title = name + ' RSS_' + rev + direc 
                            if verbose:
                                print('Extracting CSV data...')
                                print(geo_file,cal_file,dlp_file)
                            data = rss_ringoccs.ExtractCSVData(geo_file, cal_file, dlp_file)
                            #print('dir(data)',dir(data))
                            if verbose:
                                print('Finished extracting CSV data')
                            if NoOp:

                                return -1

                            if drho_km != 0:
                                data.rho_km_vals += drho_km
                            for psitype in psitypes:
                                tstart_psitype = time.time()
                                for wtype in wtypes:
                                    tstart_wtype = time.time()
                                    for res_km in res_kms:
                                        tstart_res = time.time()
                                        
                                        if verbose:
                                            print(psitype,res_km,'...')
                                        
                                        if subdivide:
                                            nranges,min_ranges,max_ranges = subdivide_inversion_range(inversion_range,drange)
                                        else:
                                            nranges = 1
                                            min_ranges = [min(inversion_range)]
                                            max_ranges = [max(inversion_range)]
                                        outfiles_to_merge = []
                                        for min_range, max_range in zip(min_ranges,max_ranges):
                                            this_inversion_range = [min_range,max_range]
                                            #print('this_inversion_range',this_inversion_range)
                                            tstart_diffrac = tstart_res 
                                            try: # may fail due to data not available for this request
                                                tau_inst = rss_ringoccs.DiffractionCorrection(
                                                   data, res_km, rng=this_inversion_range, resolution_factor=res_factor,
                                                   psitype=psitype, wtype=wtype, verbose=False)
                                                #print('dir(tau_inst)',dir(tau_inst))
                                            except:
                                                print("DiffractionCorrection failed - check data coverage, requested resolution")
                                                try: # perhaps there are some resolutions already processed, so plot them
                                                    if not plot_individual:
                                                        if xlim_run == (None,None):
                                                            xlim = inversion_range
                                                        else:
                                                            xlim = xlim_run
                                                        try: # will fail if title not defined because not a valid dlp this pass
                                                            figfile = plot_comparisons(program,results_all,plot_phase=plot_phase,y_offset=y_offset,
                                                               adjust_phase=adjust_phase,dy_offset=dy_offset,show=True,
                                                               local_path_to_data = global_path_to_data,title=title,
                                                               CORSS_8001_filepath=CORSS_8001_filepath,include_CORSS_8001=False,
                                                               figsize=figsize,fontsize_legend = fontsize_legend,xlim=xlim,ylim_tau=ylim_tau,
                                                               path_to_figs=global_path_to_local_figs,save_figfile=save_figfile)
                                                            #print('figfile = ',figfile)
                                                        except:
                                                            pass
                                                except:
                                                    pass
                                                return outfiles_complete # perhaps some of these have been successfully run
                                            tend_diffrac = time.time()
                                            if write_tau_files:
                                                rev_info = get_rev_info_from_dlp(dlp_file)
                                                #print("rev_info:",rev_info)
                                                if include_history:
                                                    tau_history = set_tau_history_from_csv(tau_inst,geo_file,cal_file,dlp_file,
                                                        tstart_diffrac,tend_diffrac,res_km,this_inversion_range,res_factor,psitype,wtype,program,
                                                        rssocc_version='1.3-beta')
                                                else:
                                                    tau_history = None
                                                if subdivide:
                                                    this_path_to_output = global_path_to_local_tmp
                                                else:
                                                    this_path_to_output = global_path_to_local_output 
                                                #print('dir(tau_inst)',dir(tau_inst))
                                                # compute tau_threshold
                                                tau_inst.tau_threshold_vals = compute_tau_threshold(cal_file,tau_inst)
                                                
                                                outfiles = write_output_files.write_output_files(tau_inst,rev_info=rev_info,
                                                        history=tau_history,local_path_to_output = this_path_to_output)
                                                outfiles_to_merge.append(outfiles)
                                                outfiles_complete.append(outfiles)
                                                if not subdivide:
                                                    add_psitype_to_taufiles(outfiles,psitype,verbose=verbose)
                                            results = {'figfile':figfile,'geo_file':geo_file,'cal_file':cal_file,'dlp_file':dlp_file,'res_km':res_km,
                                   'tau_file':None,'rkm':tau_inst.rho_km_vals,'tau':tau_from_tau_inst(tau_inst),
                                   'tau_threshold':tau_inst.tau_threshold_vals,'CORSS_8001_filepath':'None','res_km':res_km,
                                   'res_factor':res_factor,'inversion_range':this_inversion_range,'psitype':psitype,'band':band,
                                   'wtype':wtype,'Rev':rev,'direction':direc,'DSN':dsn,'kernels':'unknown'}
                                            if plot_individual:
                                                if xlim_run == (None,None):
                                                    xlim = inversion_range
                                                else:
                                                    xlim = xlim_run
                                                figfile = plot_comparisons(program,results,plot_phase=plot_phase,y_offset=y_offset,
                                                       adjust_phase=adjust_phase,dy_offset=dy_offset,show=True,
                                                       local_path_to_data = global_path_to_data,title=title,
                                                       CORSS_8001_filepath=CORSS_8001_filepath,include_CORSS_8001=False,
                                                        figsize=figsize,fontsize_legend = fontsize_legend,xlim=xlim,ylim_tau=ylim_tau,
                                                       path_to_figs=global_path_to_demo_figs,save_figfile=save_figfile)
                                                results['figfile'] = figfile
                                            else:
                                                results_all.append(results)
                                                results_complete.append(results)
                                        if subdivide and nranges > 1:
                                            outfiles_to_merge = flatten_extend(outfiles_to_merge)
                                            merge_tau_files(outfiles_to_merge,psitype,tmp2output=tmp2output,clean=clean,silent=silent)
                                        tend_res = time.time()
                                        if PRINT_TIME:
                                            print(format_time(name,rev,band,dsn,direc,dlp_res_used,psitype,res_km,tstart_res,tend_res,MINUTES=False))
                                    tend_wtype = time.time()
                                    # if PRINT_TIME:
                                    #     print('wtype '+wype+' processing time (seconds): '+f'{tend_wtype-tstart_wtype:0.2f}')
                                tend_psitype= time.time()
                                # if PRINT_TIME:
                                #     print(format_time(name,rev,band,dsn,direc,dlp_res_used,psitype,res_km,tstart_psitype,tend_psitype,MINUTES=False))
                            tend_dsn = time.time()
                            # if PRINT_TIME:
                            #     print('DSN '+dsn+' processing time (seconds): '+f'{tend_dsn-tstart_run:0.2f}')
                        tend_band= time.time()
                        # if PRINT_TIME:
                        #     print(band+'band processing time (seconds): '+f'{tend_band-tstart_dsn:0.2f}'
                    # else:
                    #     return -1 # no DLP file found
            if NoOp:
                return -1
            tend_run = time.time()
            if PRINT_TIME and outfiles_complete != []:
                print('run processing time: '+f'{(tend_run-tstart_run):0.2f} sec')
            if not plot_individual:
                #xlim = inversion_range #xlim = (np.min(tau_inst.rho_km_vals),np.max(tau_inst.rho_km_vals))
                if xlim_run == (None,None):
                    xlim = inversion_range
                else:
                    xlim = xlim_run
                try: # will fail if title not defined because not a valid dlp this pass
                    figfile = plot_comparisons(program,results_all,plot_phase=plot_phase,y_offset=y_offset,
                       adjust_phase=adjust_phase,dy_offset=dy_offset,show=True,
                       local_path_to_data = global_path_to_data,title=title,
                       CORSS_8001_filepath=CORSS_8001_filepath,include_CORSS_8001=False,
                       figsize=figsize,fontsize_legend = fontsize_legend,xlim=xlim,ylim_tau=ylim_tau,
                       path_to_figs=global_path_to_local_figs,save_figfile=save_figfile)
                    #print('figfile = ',figfile)
                except:
                    pass
        tend_feature= time.time()
        if PRINT_TIME and outfiles_complete != []:
            print(name + ' processing time: '+f'{(tend_feature-tstart_feature)/60:0.2f} min')      
    tend_all = time.time()
    if PRINT_TIME and outfiles_complete != []:
        print('Total processing time: '+f'{(tend_all-tstart_all)/60:0.2f} min')
        print(get_processor_info())
    return outfiles_complete # results_complete
    
def set_tau_history(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                    res_km,inversion_range,res_factor,psitype,wtype,program):
    """
    Construct tau_inst.history by harvesting information from geo, cal, dlp instances
    """
    set_tau_history_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                    res_km,inversion_range,res_factor,psitype,wtype,program,rssocc_version='1.3-beta')
    dlp_label_file = dlp_inst.outfiles[0]+'.LBL' # get history from DLP LBL file
    f = open(dlp_label_file, "r")
    contents = f.read().splitlines() # strips \n at end of each line
    f.close()
    for i,line in enumerate(contents):
        if line.startswith("rsr_inst history:"):
            istart = i
        if line.startswith('"'):
            istop = i
            break
    prior_history = contents[istart:istop]
    geo_file = geo_inst.outfiles[0]+'.TAB'
    cal_file = cal_inst.outfiles[0]+'.TAB'
    dlp_file = dlp_inst.outfiles[0]+'.TAB'

    try:
#        print("try setattr(tau_inst,'history',dlp_inst.history")
        setattr(tau_inst,'history',dlp_inst.history)
        tau_inst.history['Positional Args']= {'GEO file':geo_file,'CAL file':cal_file,'DLP file':dlp_file,
                                              'res_km':str(res_km)}
        tau_inst.history['Keyword Args']= {'rng':str(inversion_range),'res_factor':str(res_factor),
                'psitype':psitype,'wtype':wtype}
    
        user_name = os.getlogin()
        host_name = os.uname()[1]
        run_date = time.ctime() + ' ' + time.tzname[0]
        python_version = platform.python_version()
        operating_system = os.uname()[0]+' '+platform.platform()
    
        tau_inst.history['Run Date']=run_date
        tau_inst.history['Source File']=program
        tau_inst.history['Source Directory']=os.getcwd()
        tau_inst.history['User Name']=user_name
        tau_inst.history['Host Name']=host_name
        tau_inst.history['Python Version']=python_version
        tau_inst.history['Operating System']=operating_system
        dtminutes = (tend-tstart)/60
        sdtminutes = f'{dtminutes:0.2f} minutes'
        tau_inst.history['Additional Info']={'Prior history':prior_history,
                                             'Diffraction reconstrution time':sdtminutes}
    except:
        print('WARNING: Unable to update TAU file history')
    return tau_inst # updated tau_inst with possibly-updated history fields defined
    
def set_tau_history_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                    res_km,inversion_range,res_factor,psitype,wtype,program,rssocc_version='1.3-beta'):
    """
    Construct tau_inst.history by harvesting information from geo, cal, dlp instances
    """
    dlp_label_file = dlp_inst.outfiles[0]+'.LBL' # get history from DLP LBL file
    f = open(dlp_label_file, "r")
    contents = f.read().splitlines() # strips \n at end of each line
    f.close()
    for i,line in enumerate(contents):
        if line.startswith("rsr_inst history:"):
            istart = i
        if line.startswith('"'):
            istop = i
            break
    prior_history = contents[istart:istop]
    geo_file = geo_inst.outfiles[0]+'.TAB'
    cal_file = cal_inst.outfiles[0]+'.TAB'
    dlp_file = dlp_inst.outfiles[0]+'.TAB'

    try:
        user_name = os.getlogin()
    except:
        user_name = 'UNKNOWN'
    try:
        host_name = os.uname()[1]
    except:
        host_name = 'UNKNOWN'
    run_date = time.ctime() + ' ' + time.tzname[0]
    python_version = platform.python_version()
    operating_system = os.uname()[0]+' '+platform.platform()
    dtminutes = (tend-tstart)/60
    sdtminutes = f'{dtminutes:0.2f} minutes'
    tau_history = {'Positional Args':{'GEO file':geo_file,'CAL file':cal_file,'DLP file':dlp_file,\
        'res_km':str(res_km)},'rss_ringoccs Version': rssocc_version,\
        'Keyword Args':{'rng':str(inversion_range),'res_factor':str(res_factor),'psitype':psitype,'wtype':wtype},\
        'Run Date':run_date,'Source File':program,'Source Directory':os.getcwd(),\
        'User Name':user_name,'Host Name':host_name,'Python Version':python_version,\
        'Operating System':operating_system,\
        'Additional Info':{'Prior history':prior_history,\
        'Diffraction reconstrution time':sdtminutes}}
    return tau_history
    
def spm2date(year,doy,spmval):
    hrs = int(spmval//3600)
    mins= int(spmval//60 % 60)
    secs= int(spmval % 60)
    datestring = year + '-' + doy + 'T' + f'{hrs:02d}:{mins:02d}:{secs:02d}'
    return datestring
    
def subdivide_inversion_range(inversion_range, drange):

    min_ranges = []
    max_ranges = []

    nranges = 1
    condition = False
    while not condition:
        min_range = inversion_range[0]+(nranges-1)*drange
        min_ranges.append(min_range)
        max_range = min_range + drange
        max_range = min([max_range,inversion_range[1]])
        max_ranges.append(max_range)
        if max_range == inversion_range[1]:
            nranges -= 1
            condition=True
        nranges += 1
    return nranges,min_ranges,max_ranges
    
def tau_from_CSV(geo_file,cal_file,dlp_file,Rev,direction,DSN,res_km=1.0,res_factor=0.75,   
    inversion_range=[87475,87560],psitype='fresnel',
    wtype='kbmd20',CORSS_8001_filepath=None,include_CORSS_8001=False,plot=True,local_path_to_data =global_path_to_data,
    y_offset=0,dy_offset=0, local_path_to_demo_figs=global_path_to_demo_figs,save_figfile=True,
    figsize=(8,5),fontsize_legend = 'x-small',
    plot_phase = False, adjust_phase = False,title=None,
    verbose=True,show=True,program='demo_tau_from_CSV'):

    is_res_km_valid(dlp_file,res_factor,res_km) # returns exception if invalid requesed res_km
    data = rss.ExtractCSVData(geo_file, cal_file, dlp_file)
    tau_inst = rss.DiffractionCorrection(data, res_km, rng=inversion_range, res_factor=res_factor,
           psitype=psitype, wtype=wtype, verbose=verbose)
#added 2026 Jan 02
    tau_inst.tau_threshold_vals = compute_tau_threshold(cal_file,tau_inst)
    #print(dir(tau_inst))
    figfile = None
    results = {'figfile':figfile,'geo_file':geo_file,'cal_file':cal_file,'dlp_file':dlp_file,
               'res_km':res_km,
               'tau_file':None,'rkm':tau_inst.rho_km_vals,'tau':tau_from_tau_inst(tau_inst),
               'CORSS_8001_filepath':CORSS_8001_filepath,'res_km':res_km,
               'res_factor':res_factor,'inversion_range':inversion_range,'psitype':psitype,
               'wtype':wtype,'Rev':Rev,'direction':direction,'DSN':DSN,'kernels':kernels}
    if plot:
        figfile = plot_comparisons(program,results,plot_phase=plot_phase,y_offset=y_offset,
                                   adjust_phase=adjust_phase,dy_offset=dy_offset,show=show,
                                   local_path_to_data = local_path_to_data,title=title,
                                   CORSS_8001_filepath=CORSS_8001_filepath,include_CORSS_8001=include_CORSS_8001,
                                   figsize=figsize,fontsize_legend = fontsize_legend,
                                   path_to_figs=global_path_to_demo_figs,save_figfile=save_figfile)
        results['figfile'] = figfile
    return results
    
def tau_from_tau_inst(tau_inst):
    '''
    Compute optical depth from tau_inst
    '''
    pwr = np.abs(tau_inst.T_out)**2
    mu = np.sin(np.radians(np.abs(tau_inst.B_deg_vals)))
    tau = -mu * np.log(pwr)
    return tau 
    
# integrand of tau
def tau_int(a,Q,q=3.1):
    return a**2. * Q * n(a,q)
    
def trim_dlp_file(dlp_file,processing_range,path_to_output=global_path_to_local_output,overwrite=False,write=True,verbose=False):
    dirname = os.path.dirname(dlp_file)
    basename = os.path.basename(dlp_file)
    output_dir = path_to_output + dirname[dirname.index('Rev'):]+'/'
    os.makedirs(output_dir,exist_ok=True)
    insert_string = str(int(processing_range[0]))+'-'+str(int(processing_range[1]))
    trimmed_dlp_file = output_dir + basename.replace('M_','M_'+insert_string+'_')
    path = trimmed_dlp_file
    exists = os.path.exists(trimmed_dlp_file)
    # if verbose and exists:
    #     print('Trimmed DLP already exists')
#Don't write the file if it already exists unless overwrite=True
    if write and (not os.path.exists(path) or overwrite):
        rkm = np.loadtxt(dlp_file,delimiter = ',',usecols=(0))
        L = np.where((rkm >= processing_range[0]) & (rkm <= processing_range[1]))[0]
        with open(dlp_file) as f:
            contents = f.readlines()
        file = open(path,'w')
        for LL in L:
            line = contents[LL]
            file.write(line)
        file.close()
        if verbose:
            print('wrote '+path)
        orig_dlp_label_file = dlp_file.replace('TAB','LBL')
        new_dlp_label_file = path.replace('TAB','LBL')
        shutil.copy(orig_dlp_label_file,new_dlp_label_file)
        if verbose:
            print('copied',os.path.basename(orig_dlp_label_file),' to ',new_dlp_label_file)
    # else:
    #     if verbose:
    #         print('write and (not os.path.exists(path) or overwrite)',write and (not os.path.exists(path) or overwrite))
    return trimmed_dlp_file
        
def update_file_version(pathtofile):
#    print('updating version number')
#    print('input filename',pathtofile)
    newpathtofile = pathtofile
    lenf = len(pathtofile)
    index=lenf-9
    version = int(pathtofile[lenf-8:lenf-4])
    print('version:',version)
#    print('os.path.exists(newpathtofile)',os.path.exists(newpathtofile))
    while os.path.exists(newpathtofile):
        version += 1
        sversion = f'{version:04d}'
        newpathtofile = newpathtofile[0:index+1]+sversion+newpathtofile[index+5:]
#        print(newpathtofile)
#        print(os.path.exists(newpathtofile))
    return newpathtofile

def update_taufile_radius_corrections(taufile,update_dr_pole=True,update_dr_trajectory=True,
                                     fit_number=1,NMAX_POLE = 500,NMAX_TRAJECTORY = 1000,
                                     write_only_if_zero=True,verbose=False,
                                     save_original=True):
    with open(taufile,'r') as f:
        contents = f.readlines()
    if write_only_if_zero:
        dr_pole_,dr_traj_=np.loadtxt(taufile,delimiter=',',usecols=[1,2]).T
        if np.any(dr_pole_) or np.any(dr_traj_):
            if verbose:
                print('taufile not overwritten! dr_pole,dr_traj already non-zero')
                print('No need to copy original to original~')
            return
    if save_original:
        shutil.copyfile(taufile, taufile+'~')
        if verbose:
            print('Copied original to original~')
            print(taufile+'~')
            print('Likely moved to tmp*/old subdirectory')
    if update_dr_pole:
        dr_pole,void,void = radius_correction_pole(taufile,fit_number=fit_number,NMAX=NMAX_POLE)
        assert len(contents) == len(dr_pole),"len(contents) != len(dr_pole)"
    if update_dr_trajectory:
        dr_trajectory,void,void = radius_correction_trajectory(taufile,NMAX=NMAX_TRAJECTORY)
        assert len(contents) == len(dr_pole),"len(contents) != len(dr_pole)"
    ofile = taufile
    with open (ofile,'w',newline='\r\n') as f:
        for i,line_ in enumerate(contents):
            if verbose and i<5:
                print(line_)
            line = list(line_)
            if update_dr_pole:
    #                '  87482.500000,  0.152845,  0.941343, '
    #            print(line)
                line[16:25]=f'{dr_pole[i]:9.6f}'
    #            print(line)
            if update_dr_trajectory:
    #                '  87482.500000,  0.152845,  0.941343, '
    #            print(line)
                line[27:36]=str(f'{dr_trajectory[i]:9.6f}')
    #            print(line)
            newline = "".join(line)
            if verbose and i<5:
                print(newline)

            # assert i<5,'Halt after 5'
                
            f.write(newline)
               
# write TAU *.LBL and *.TAB files, including version with psitype appended to root name
def write_tau_files(tau_inst,geo_inst,cal_inst,dlp_inst,dlp_file,tstart,tend,
                   res_km,inversion_range,res_factor,psitype,wtype,program,
                   local_path_to_output='../output/',verbose=True):
    """
    write TAU *.LBL and *.TAB files, including version with psitype appended to root name
    tau_inst,geo_inst,cal_inst,dlp_inst: required instances of tau, geo, cal, dlp
    tstart,tend: SPM of start and end times
    res_km: recontruction resolution
    inversion_range: (rmin, rmax) radial range in km of reconstruction range of inversion
    res_factor: reconstruction scale factor to match PDS convention
    wtype: processing window type
    program: program name
    local_path_to_output: local path to rss_ringoccs/output/ directory

    returns the path to the TAU TAB file created.
    """
# update the history in the label file if possible
    
#    print('write_tau_files:local_path_to_output',local_path_to_output)

    set_tau_history(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                    res_km,inversion_range,res_factor,psitype,wtype,program)
# write the file - need to add keyword rev_info since missing from 
    rev_info = dlp_inst.rev_info
    outfiles = rss.tools.write_output_files.write_output_files(tau_inst,rev_info=rev_info,
            local_path_to_output=local_path_to_output)

# copy output files with psitype added
    for suffix in ['.LBL','.TAB']:
        sourcefile = outfiles[0]+suffix
        destfile = outfiles[0]+psitype+suffix
        shutil.copy(sourcefile,destfile)
        if verbose:
            print('Copied\n',os.path.basename(sourcefile),'\nto\n',os.path.basename(destfile))
    return destfile # return the tau_file 

def write_tau_files_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                   res_km,inversion_range,res_factor,psitype,wtype,program,
                   local_path_to_output=global_path_to_output,verbose=True):
    """
    write TAU *.LBL and *.TAB files, including version with psitype appended to root name
    tau_inst,geo_inst,cal_inst,dlp_inst: required instances of tau, geo, cal, dlp
    tstart,tend: SPM of start and end times
    res_km: recontruction resolution
    inversion_range: (rmin, rmax) radial range in km of reconstruction range of inversion
    res_factor: reconstruction scale factor to match PDS convention
    wtype: processing window type
    program: program name
    local_path_to_output: local path to rss_ringoccs/output/ directory

    returns the path to the TAU TAB file created.
    """

    tau_history=set_tau_history_from_inst(tau_inst,geo_inst,cal_inst,dlp_inst,tstart,tend,
                    res_km,inversion_range,res_factor,psitype,wtype,program)
#    print('prior to call of write_output_files: tau_history',tau_history)
# write the file - need to add keyword rev_info since missing from tau_inst.
    rev_info = dlp_inst.rev_info
    outfiles = rss.tools.write_output_files.write_output_files(tau_inst,rev_info=rev_info,
            history=tau_history,local_path_to_output=local_path_to_output)

# copy output files with psitype added
    for suffix in ['.LBL','.TAB']:
        sourcefile = outfiles[0]+suffix
        destfile = outfiles[0]+'_'+psitype+suffix
        shutil.copy(sourcefile,destfile)
        if verbose:
            print('Copied\n',os.path.basename(sourcefile),'\nto\n',os.path.basename(destfile))
    return destfile # return the tau_file
