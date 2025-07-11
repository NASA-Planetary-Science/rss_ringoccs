"""

write_uranus_output_files.py

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
import os
from time import strftime
import rss_ringoccs as rss
from .pds3_uranus_geo_series import write_uranus_geo_series
from .pds3_tau_series import write_tau_series

chord_revnums = ['053', '054', '056', '057', '058', '060', '063', '064',
                    '067', '079', '081', '082', '084', '089', '253']

func_typ = {'GEO': write_uranus_geo_series, 'TAU': write_tau_series}

def write_uranus_output_files(inst):
    """
    Write output (geo, cal, dlp, tau) *.TAB and *.LBL files, depending on
    instance given.

    Args:
        inst (instance):
            Instance of either Geometry, Calibration, NormDiff, or
            DiffractionCorrection classes.
    """
    rev_info = inst.rev_info

    # Check for instance type and write that specific file
    if isinstance(inst, rss.occgeo.occgeo_uranus.Uranus_Geometry):
        filtyp = 'GEO'

    elif isinstance(inst, rss.DiffractionCorrection):
        filtyp = 'TAU_' + str(int(inst.input_resolution_km * 1000)).zfill(5) + 'M'
    else:
        raise TypeError('invalid instance for write_output_files.')

    construct_output_filename(rev_info, inst, filtyp)
    return None

def construct_filepath(rev_info, filtyp):
    title_out = []
    outdir_out = []

    pd1 = rev_info['prof_dir']
    pd1 = [pd1]
    pd2 = ''

    doy = rev_info['doy']
    band = rev_info['band']
    year = rev_info['year']
    dsn = rev_info['dsn'].split('-')[-1]
    rev = rev_info['rev_num']
    ring = rev_info['ring']

    if rev in chord_revnums:
        if 'DLP' in filtyp or 'TAU' in filtyp or 'Summary' in filtyp:
            pd2 = 'C'

    for dd in pd1:
        filestr = ('VGR2_' + str(band) + str(dsn) + '_' + dd + '_URING_'
                    + ring )

        dirstr = ('../output/' + dd + '/' + ring + '/')

        # Create output file name without file extension
        curday = strftime('%Y%m%d')
        out1 = dirstr + filestr + '_' + filtyp.upper() + '_' + curday

        if os.path.isdir(dirstr):

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
            print('\tCreating directory:\n\t\t' + dirstr)
            #os.mkdir(dirstr)
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
    for n in range(ndirs):
        title = titles[n]
        outdir = outdirs[n]
        func_typ[filtyp[0:3]](rev_info, inst, title, outdir,
                rev_info['prof_dir'])

    return None

"""
Revisions:
"""
