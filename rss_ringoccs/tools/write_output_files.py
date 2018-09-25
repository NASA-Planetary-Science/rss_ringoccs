"""

write_output_files.py

Purpose: Write output *.TAB data and corresponding *.LBL label file.

Revisions:
    2018 Sep 19 - jfong - original
    2018 Sep 24 - jfong - hardcode relative path (../output/*)
"""
import sys
from .pds3_geo_series import write_geo_series
from .pds3_cal_series import write_cal_series
from .pds3_dlp_series import write_dlp_series
from .pds3_tau_series import write_tau_series
from time import strftime

sys.path.append('../../')
import rss_ringoccs as rss
sys.path.remove('../../')

import pdb
import os
func_typ = {'GEO': write_geo_series,
        'CAL': write_cal_series,
        'DLP': write_dlp_series,
        'TAU': write_tau_series}

def write_output_files(inst):
    rev_info = inst.rev_info
    
    # Check for instance type and write that specific file
    if isinstance(inst, rss.occgeo.Geometry):
        filtyp = 'GEO'

    elif isinstance(inst, rss.calibration.Calibration):
        filtyp = 'CAL'

    elif isinstance(inst, rss.calibration.NormDiff):
        filtyp = 'DLP_' + str(int(inst.dr_km * 1000 * 2)).zfill(4) + 'M'

    elif isinstance(inst, rss.diffcorr.DiffractionCorrection):
        filtyp = 'TAU_' + str(int(inst.res * 1000)).zfill(5) + 'M'
    else:
        print('invalid instance!')

    construct_output_filename(rev_info, inst, filtyp)
    return None

def construct_output_filename(rev_info, inst, filtyp):
    pd1 = (rev_info['prof_dir'].split('"')[1])[0]
    if pd1 == 'B':
        pd1 = ['I','E']
        pd2 = 'C'
    else:
        pd1 = [pd1]
        pd2 = ''

    doy = rev_info['doy']
    band = rev_info['band'].split('"')[1]
    year = rev_info['year']
    dsn = rev_info['dsn'].split('-')[-1]
    rev = rev_info['rev_num']

    for dd in pd1:
        filestr = ('RSS_' + str(year) + '_' + str(doy) + '_' + str(band) +
            str(dsn) + '_' + dd)

        dirstr = ('../output/Rev' + rev + '/' + dd + '/' + 'Rev' + rev +
                pd2 + dd + '_' + filestr + '/')

#    os.system('[ ! -d ' + output_directory + ' ] && mkdir -p '
#            + output_directory)
#    fsearch = dirstr + filestr + '*' + filtyp.upper() + '*' + '.*'

        #cwds = (os.getcwd()).split('/')
        # index of first appearance of rss_ringoccs
        #ind = cwds.index('rss_ringoccs')
        #relpath = '../' * int(len(cwds)-ind-1)

        # Create output file name without file extension
        out1 = dirstr + filestr + '_' + filtyp.upper() + strftime('_%Y%m%d')

        if os.path.exists(dirstr):

            # Check for most recent file and order them by date
            dirfiles = [x for x in os.listdir(dirstr) if
                    x.startswith('.')==False]


            if len(dirfiles) == 0:
                seq_num = '0001'
            else:
                sfn = [x.split('_') for x in dirfiles]
                sqn0 = [(x[-2]+x[-1][0:4]) for x in sfn if (x[-3]==filtyp)]
                if len(sqn0) == 0:
                    seq_num = '0001'
                else:
                    sqn1 = [int(x) for x in sqn0]

                    seq_num = (str(sqn1[-1] + 1))[-4:]

            out2 = out1 + '_' + seq_num



        else:
            print('\tCreating directory:\n\t\t' + dirstr)
            os.system('[ ! -d ' + dirstr + ' ] && mkdir -p ' + dirstr)

            seq_num = '0001'
            out2 = out1 + '_' + seq_num

        title = out2.split('/')[-1]
        outdir = '/'.join(out2.split('/')[0:-1]) + '/'

        func_typ[filtyp[0:3]](rev_info, inst, title, outdir, 
                rev_info['prof_dir'])

        
    return None
