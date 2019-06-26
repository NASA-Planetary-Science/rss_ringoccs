"""

write_output_files.py

:Purpose: 
    Functions relating to writing an output file.

:Dependencies:
    #. sys
    #. pds3_geo_series
    #. pds3_cal_series
    #. pds3_dlp_series
    #. pds3_tau_series
    #. time
    #. os

"""
import sys
from .pds3_geo_series import write_geo_series
from .pds3_cal_series import write_cal_series
from .pds3_dlp_series import write_dlp_series
from .pds3_tau_series import write_tau_series
#from .pds3_spectro_image import write_spectro_image
import time
from time import strftime

sys.path.append('../../')
import rss_ringoccs as rss
sys.path.remove('../../')

import os


func_typ = {'GEO': write_geo_series,
        'CAL': write_cal_series,
        'DLP': write_dlp_series,
        'TAU': write_tau_series}
        #'SPECTRO': write_spectro_image}

def write_output_files(inst, add_text=None):
    """
    Write output (geo, cal, dlp, tau) *.TAB and *.LBL files, depending on 
    instance given.

    Args:
        inst (instance):
            Instance of either Geometry, Calibration, NormDiff, or
            DiffractionCorrection classes.
    """
    if add_text is not None:
        add = add_text
    else:
        add = ''
    rev_info = inst.rev_info

    # Check for instance type and write that specific file
    if isinstance(inst, rss.occgeo.Geometry):
        filtyp = 'GEO' + add

    elif isinstance(inst, rss.calibration.Calibration):
        filtyp = 'CAL' + add

    elif isinstance(inst, rss.calibration.DiffractionLimitedProfile):
        filtyp = 'DLP_' + str(int(inst.dr_km * 1000 * 2)).zfill(4) + 'M' + add

    elif isinstance(inst, rss.diffrec.DiffractionCorrection):
        filtyp = 'TAU_' + str(int(inst.input_res * 1000)).zfill(5) + 'M' + add
    #elif isinstance(inst, rss.scatter.Scatter):
    #    filtyp = 'SCATTER_' + inst.band + (inst.dsn).split('-')[1] + add
    else:
        print('invalid instance!')

    outfiles = construct_output_filename(rev_info, inst, filtyp)
    return outfiles

def construct_filepath(rev_info, filtyp):
    title_out = []
    outdir_out = []

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


    if 'DIR' in rev_info:
        if 'DLP' in filtyp or 'TAU' in filtyp or 'Summary' in filtyp:
            pd2 = 'C'

    if 'PER' in rev_info:
        pd2 = 'P' + pd2

    for dd in pd1:
        filestr = ('RSS_' + str(year) + '_' + str(doy) + '_' + str(band) +
            str(dsn) + '_' + dd)

        dirstr = ('../output/Rev' + rev + '/Rev' + rev + pd2 + dd 
                + '/' + 'Rev' + rev + pd2 + dd + '_' + filestr + '/')

        # Create output file name without file extension
        curday = strftime('%Y%m%d')
        out1 = dirstr + filestr + '_' + filtyp.upper() + '_' + curday

        if os.path.exists(dirstr):

            # Check for most recent file and order them by date
            dirfiles = [x for x in os.listdir(dirstr) if
                    x.startswith('.')==False]


            if len(dirfiles) == 0:
                seq_num = '0001'
            else:
                sqn0 = [x.split('_')[-2]+x.split('_')[-1][0:4] for x in dirfiles
                        if (filtyp.upper() in x) and (x.split('_')[-2]==curday)]
                if len(sqn0) == 0:
                    seq_num = '0001'
                else:
                    sqn1 = [int(x) for x in sqn0]

                    seq_num = (str(sqn1[-1] + 1))[-4:]

            out2 = out1 + '_' + seq_num



        else:
            os.system('[ ! -d ' + dirstr + ' ] && mkdir -p ' + dirstr)

            seq_num = '0001'
            out2 = out1 + '_' + seq_num

        title = out2.split('/')[-1]
        outdir = '/'.join(out2.split('/')[0:-1]) + '/'
        title_out.append(title)
        outdir_out.append(outdir)
    return title_out, outdir_out
        

def construct_output_filename(rev_info, inst, filtyp):

    titles, outdirs = construct_filepath(rev_info, filtyp)
    ndirs = len(titles)
    outfiles = []
    for n in range(ndirs):
        title = titles[n]
        outdir = outdirs[n]
        outfiles.append(outdir + title)
        func_typ[filtyp[0:3]](rev_info, inst, title, outdir,
                rev_info['prof_dir'])

    return outfiles

"""
Revisions:
    2019 Apr 09 - jfong - add Scatter outfile
"""
