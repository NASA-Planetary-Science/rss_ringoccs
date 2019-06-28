### BEGIN USER INPUT
# raw data filepath
rawpath = '/Volumes/sflury001/Research/TC2017/data/VoyagerUranus/data/raw/'#'../data/'
# global setting of verbose mode
verbose = True              # optional verbose mode
# global kernel file
kernels = '../tables/VG2_URING_kernel.ker' # kernel file
# reference files for DLP processing
geo_date = '20190628'        # GEO file date
geo_sn = '0001'              # GEO file serial number
cal_date = '20190628'        # CAL file date
cal_sn = '0001'              # CAL file serial number
# reference info for DLP processing
dr_km_desired = 0.005       # DLP sampling resolution
# reference files for Reconstruction
dlp_res = '005'             # DLP resolution in meters
dlp_date = '20190522'       # DLP file date
dlp_sn = '0002'             # DLP file serial number
# settings for Reconstruction processing
window = 'kb25'             # reconstruction window for Fresnel transform
res = 0.02                  # desired reconstruction resolution
psitype = 'cfresnel4'       # method for approximating psi
### END USER INPUT
#
### ANCILARY INFORMATION
# ring names
rnames = ['6','5','4','A','B','N','G','D','E']
# dictionary of ring locations in SPM
rings_spm = {
    'I':{'6':82001.23199844, '5':81947.86100006, '4':81916.74900055,
         'A':81648.5359993, 'B':81537.59500122, 'N':81349.53099823,
         'G':81295.03699875, 'D':81211.89100075, 'L':81001.20190001,
         'E':80906.92000008},
    'E':{'6':4703.43400002,'5':4746.70200014, '4':4795.93099976,
         'A':5055.85999966, 'B':5175.7670002 , 'N':5357.72299957,
         'G':5412.36200047,'D':5493.19499969, 'L':5701.30999994,
         'E':5859.81299973}}
# dictionary of ring locations in KM
rings_km = {
    'E':{'6':41794.93, '5':42155.17, '4':42555.67,
         'A':44686.59, 'B':45673.33, 'N':47176.51,
         'G':47628.26, 'D':48297.35, 'E':51342.23},
    'I':{'6':41871.35, '5':42304.04, '4':42556.28,
         'A':44736.75, 'B':45640.74, 'N':47176.28,
         'G':47621.59, 'D':48300.57, 'E':50796.85}}
# dictionary of ring widths in KM
ring_widths = {
    'I':{'6':1.52, '5':2.75, '4':1.95,'A':10.59, 'B':7.03, 'N':1.54,
         'G':3.83, 'D':6.7, 'E':22.43},
    'E':{'6':1.72, '5':2.62, '4':2.67, 'A':4.22, 'B':11.19, 'N':1.53,
         'G':1.63, 'D':2.7, 'E':74.93}}
