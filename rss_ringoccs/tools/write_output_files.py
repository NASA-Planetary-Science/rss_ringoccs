
"""
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
import rss_ringoccs as rss
import os


func_typ = {'GEO': write_geo_series,
        'CAL': write_cal_series,
        'DLP': write_dlp_series,
        'TAU': write_tau_series}
        #'SPECTRO': write_spectro_image}

def write_output_files(inst, add_text=None,rev_info=None,local_path_to_tables='../tables/',
                       add_suffix=None,
                       history=None,local_path_to_output = '../output/',verbose=False):
    """
    Write output (geo, cal, dlp, tau) *.TAB and *.LBL files, depending on 
    instance given.

    Arguments
        :inst (*object*): Instance of either Geometry, Calibration, 
                         DiffractionLimitedProfile, or
                         DiffractionCorrection

    Keyword Arguments
        :add_text (*str*): Additional string to be added to filename.
        :rev_info (*object*): use this rev_info instead of inst.rev_info
        :history (*object*): use this history instead of inst.history
        
    """
    if verbose:
        print('write_output_files: local_path_to_output:',local_path_to_output)
        print('write_output_files: local_path_to_tables:',local_path_to_tables)
    if add_text is not None:
        add = add_text
    else:
        add = ''
    if add_suffix is not None:
        add_sfx = add_suffix
    else:
        add_sfx = None

    # Check for instance type and write that specific file
    if isinstance(inst, rss.occgeo.Geometry):
        filtyp = 'GEO' + add

    elif isinstance(inst, rss.calibration.Calibration):
        filtyp = 'CAL' + add

    elif isinstance(inst, rss.calibration.DiffractionLimitedProfile):
        filtyp = 'DLP_' + str(int(inst.dr_km * 1000 * 2)).zfill(4) + 'M' + add

    elif isinstance(inst, rss.DiffractionCorrection):
        filtyp = 'TAU_' + str(int(inst.input_res * 1000)).zfill(5) + 'M' + add
        if verbose:
            print('TAU filtyp:',filtyp)
    #elif isinstance(inst, rss.scatter.Scatter):
    #    filtyp = 'SCATTER_' + inst.band + (inst.dsn).split('-')[1] + add
    else:
        raise TypeError('invalid instance for write_output_files.')

    if rev_info == None: # use inst.rev_info unless supplied by call
        rev_info = inst.rev_info
    if verbose:
        print('rev_info:',rev_info)
    outfiles = construct_output_filename(rev_info, inst, filtyp,history=history,local_path_to_output=local_path_to_output,add_suffix=add_sfx,verbose=verbose)
    if verbose:
        print('outfiles:',outfiles)
        print('add:',add)
    return outfiles

def construct_filepath(rev_info, filtyp,local_path_to_output = '../output/',add_suffix=None, verbose=False):
    """
    Construct output filepath 

    Arguments
        :rev_info (*dict*): Dictionary with information related to ring
                            occultation.
        :filtyp (*str*): Type of file (e.g., 'GEO', 'CAL, 'DLP', 'TAU')

    Returns
        :title_out (*list*): List of strings of output file names
        :outdir_out (*list*): List of strings of output directory path
    """
    title_out = []
    outdir_out = []

# this does not work universally. If profile direction is BOTH and occ direction
# is INGRESS or EGRESS, this would add a C for Chord to the tau file
# need to know if this is a chord occultation or not
# Hard-wiring this based on rev number
    #verbose = True # temporary
    #print('rev_info:',rev_info)
    if 'TAU' in filtyp:
        pd1 = (rev_info['prof_dir'].split('"')[1])[0]
        rev_number_int = int(rev_info['rev_num'])
        if (rev_number_int >= 53) and (rev_number_int <= 89):
            pd2 = 'C'
        else:
            pd2 = ''
    else:
        pd1 = (rev_info['prof_dir'].split('"')[1])[0]
        #pd1 = (rev_info['occ_dir'].split('"')[1])[0]

        if pd1 == 'B':
            pd1 = ['I','E']
            pd2 = 'C'
        else:
            pd1 = [pd1]
            pd2 = ''
    if verbose:
        print('construct_filepath: pd1,pd2:',pd1,pd2)


    doy = rev_info['doy']
    band = rev_info['band'].split('"')[1]
    year = rev_info['year']
    if int(year)>2010:
        dsn = rev_info['dsn'].split('-')[-1] + rev_info['ul_dsn'].split('-')[-1]
    else:
        dsn = rev_info['dsn'].split('-')[-1]
    rev = rev_info['rev_num']
    if verbose:
        print('construct_filepath: dsn,rev,local_path_to_output',dsn,rev,local_path_to_output)

# 2025 Apr 23 - the following two entries are not in the standard Cassini rev_info
# These may be Voyager Canberra and Perth

    if 'DIR' in rev_info:
        if 'DLP' in filtyp or 'TAU' in filtyp or 'Summary' in filtyp:
            pd2 = 'C'

    if 'PER' in rev_info:
        pd2 = 'P' + pd2
    if verbose:
        print('construct_filepath: now pd2 =',pd2)

    for dd in pd1:
        filestr = ('RSS_' + str(year) + '_' + str(doy) + '_' + str(band) +
            str(dsn) + '_' + dd)

        dirstr = (local_path_to_output+'Rev' + rev + '/Rev' + rev + pd2 + dd 
                + '/' + 'Rev' + rev + pd2 + dd + '_' + filestr + '/')
        if verbose:
            print('construct_filepath: dd,filestr,dirstr =',dd,filestr,dirstr)

        # Create output file name without file extension
        curday = strftime('%Y%m%d')
        out1 = dirstr + filestr + '_' + filtyp.upper() + '_' + curday
        if verbose:
            print('construct_filepath: curday,out1 =',curday,out1)

        if os.path.exists(dirstr):

            # Check for most recent file and order them by date
            dirfiles = [x for x in os.listdir(dirstr) if
                    x.startswith('.')==False]
            if verbose:
                print('construct_filepath: path exists,len(dirfiles, dirfiles = ',len(dirfiles),dirfiles)


            if len(dirfiles) == 0:
                seq_num = '0001'
            else:
                sqn0 = [x.split('_')[-2]+x.split('_')[-1][0:4] for x in dirfiles
                        if (filtyp.upper() in x) and (x.split('_')[-2]==curday)]
                if verbose:
                    print('construct_filepath: sqn0:',sqn0)
                if len(sqn0) == 0:
                    seq_num = '0001'
                else:
                    sqn1 = [int(x) for x in sqn0]
                    sqn1max = max(sqn1)
# BUG sqn1 is not sorted!
#                    seq_num = (str(sqn1[-1] + 1))[-4:]
                    seq_num = (str(sqn1max + 1))[-4:]
                    if verbose:
                        print('construct_filepath: sqn1, sqn1max,seq_num:',sqn1,type(sqn1),sqn1max,seq_num)

            out2 = out1 + '_' + seq_num
            if verbose:
                print('construct_filepath: out2,out1,seq_num:',out2,out1,seq_num)

        else:
            os.system('[ ! -d ' + dirstr + ' ] && mkdir -p ' + dirstr)
            if verbose:
                print('construct_filepath: mkdir -p dirstr:',dirstr)

            seq_num = '0001'
            out2 = out1 + '_' + seq_num
            if verbose:
                print('construct_filepath: out2, out1, seq_num:',out2,out1,seq_num)

        title = out2.split('/')[-1]
        outdir = '/'.join(out2.split('/')[0:-1]) + '/'
        title_out.append(title)
        outdir_out.append(outdir)
        if verbose:
            print('construct_filepath: title,outdir,title_out,outdir_out\n',title,outdir,title_out,outdir_out)
    return title_out, outdir_out
        

def construct_output_filename(rev_info, inst, filtyp,history=None,add_suffix=None,local_path_to_output='../output/',verbose=False):
    """
    Construct output filepath 

    Arguments
        :rev_info (*dict*): Dictionary with information related to ring
                            occultation.
        :inst (*object*): Instance of either Geometry, Calibration, 
                         DiffractionLimitedProfile, or
                         DiffractionCorrection
        :filtyp (*str*): Type of file (e.g., 'GEO', 'CAL, 'DLP', 'TAU')

    Returns
        :outfiles (*list*): List of output filepaths
    """

    titles, outdirs = construct_filepath(rev_info, filtyp,local_path_to_output=local_path_to_output,verbose=verbose)
    if verbose:
        print('construct_output_filename: titles, outdirs:',titles,outdirs)
    ndirs = len(titles)
    outfiles = []
    for n in range(ndirs):
        title = titles[n]
        outdir = outdirs[n]
        outfiles.append(outdir + title)
        filetype = filtyp[0:3]
        if filetype != 'TAU':
            func_typ[filetype](rev_info, inst, title, outdir,rev_info['prof_dir'])
        else:
            #print('write_output_files: history=',history)
#            print('TP1:add_suffix',add_suffix)
            func_typ['TAU'](rev_info, inst, title, outdir,rev_info['prof_dir'],history=history,add_suffix=add_suffix)

    return outfiles

"""
Revisions:
"""
