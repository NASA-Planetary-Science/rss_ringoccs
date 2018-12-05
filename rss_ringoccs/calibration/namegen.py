#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from time import strftime

sys.path.append('../..')
import rss_ringoccs as rss
sys.path.remove('../..')

# make output TAB/LBL file name
def tabname(rev_info,filetyp):
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
                sfn = [x.split('_') for x in dirfiles]
                sqn0 = [(x[-2]+x[-1][0:4]) for x in sfn if (x[-3]==filtyp)
                        and (x[-2]==curday)]
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

    return title,outdir

# make output plot file name following *.TAB and *.LBL files
def plotname(rev_info,filetype):
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

        # Create output file name without file extension
        curday = strftime('%Y%m%d')
        out1 = dirstr + filestr + '_' + filetype + '_' + curday

        if os.path.exists(dirstr):

            # Check for most recent file and order them by date
            dirfiles = [x for x in os.listdir(dirstr) if
                    x.startswith('.')==False]


            if len(dirfiles) == 0:
                seq_num = '0001'
            else:
                sfn = [x.split('_') for x in dirfiles]
                sqn0 = [(x[-2]+x[-1][0:4]) for x in sfn if (x[-3]==filetype)
                        and (x[-2]==curday)]
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

        return out2+'.PNG'

# make output plot file name following *.TAB and *.LBL files
def CSVname(rev_info,filetype):
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

        # Create output file name without file extension
        curday = strftime('%Y%m%d')
        out1 = dirstr + filestr + '_' + filetype + '_' + curday

        if os.path.exists(dirstr):

            # Check for most recent file and order them by date
            dirfiles = [x for x in os.listdir(dirstr) if
                    x.startswith('.')==False]


            if len(dirfiles) == 0:
                seq_num = '0001'
            else:
                sfn = [x.split('_') for x in dirfiles]
                sqn0 = [(x[-2]+x[-1][0:4]) for x in sfn if (x[-3]==filetype)
                        and (x[-2]==curday)]
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

        return out2+'.CSV'
