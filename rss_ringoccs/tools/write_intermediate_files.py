"""

write_intermediate_files.py

Purpose: This function writes output intermediate data files (such as
         frequency offset text file, frequency residual fit pickle file,
         power normalization fit pickle file) to a specific
         directory. If the specified intermediate file already exists, a new
         intermediate file with the current date will be added.
         If the output directory does not exist, it will also
         be created. Please refer to pXX of the User's Guide for our
         directory structure nomenclature.

Revisions:
    2018 Sep 18 - jfong - original
    2018 Sep 19 - jfong - add 4-digit sequence number after date to outfile str
    2018 Sep 20 - jfong - search if file type exists for seq_num='0001'
"""

import os
from time import strftime
import pdb
import pickle
from .date_to_rev import date_to_rev
import numpy as np


def write_intermediate_files(year, doy, band, dsn, profdir, filtyp,
        contents):
    # Remove 'DSS-' from dsn name
    dsn = dsn.split('-')[-1]

    # Retrieve rev number in 3-digit format
    rev = date_to_rev(year, doy)

    # Rename profdir to be I for ingress, E for egress, CI/CE for chord
    pd1 = (profdir.split('"')[1])[0]
    if pd1 == 'B':
        pd1 = ['I','E']
        pd2 = 'C'
    else:
        pd1 = [pd1]
        pd2 = ''


    # Construct directory structure nomenclature
    #   RevXXX/P/RevXXXP_RSS_YEAR_DOY_BDD_P/RSS_YEAR_DOY_BDD_filtype_DATE_#
    # sample chord: Rev054/E/Rev054CE_RSS_2007_353_K55_E/
    # Rev054/E/Rev054CE_RSS_2007_353_K55_E
#    pdb.set_trace()

    for dd in pd1:
        filestr = ('RSS_' + str(year) + '_' + str(doy) + '_' + str(band) +
            str(dsn) + '_' + dd)
    
        dirstr = ('output/Rev' + rev + '/' + dd + '/' + 'Rev' + rev + 
                pd2 + dd + '_' + filestr + '/')

#    os.system('[ ! -d ' + output_directory + ' ] && mkdir -p '
#            + output_directory)
#    fsearch = dirstr + filestr + '*' + filtyp.upper() + '*' + '.*'

        cwds = (os.getcwd()).split('/')
        # index of first appearance of rss_ringoccs
        ind = cwds.index('rss_ringoccs')
        relpath = '../' * int(len(cwds)-ind-1)

        # Create output file name without file extension
        dirsrch = relpath + dirstr
        out1 = dirsrch + filestr + '_' + filtyp.upper() + strftime('_%Y%m%d')

        # Check if directory exists, if not, create it
        if os.path.exists(dirsrch):

            # Check for most recent file and order them by date
            dirfiles = [x in os.listdir(dirsrch) if not x.startswith('.')]
            
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
            print('\tCreating directory:\n\t\t' + dirsrch)
            os.system('[ ! -d ' + dirsrch + ' ] && mkdir -p ' + dirsrch)

            seq_num = '0001'
            out2 = out1 + '_' + seq_num

        if filtyp == 'FOF':
            outfile = out2 + '.TXT'
            f_spm = contents['f_spm']
            f_offset = contents['f_offset']
            print('\tSaving frequency offset calculations to:\n\t\t'
                    + outfile)
            np.savetxt(outfile, np.c_[f_spm, f_offset], fmt='%32.16F'*2)
            
            # check data contents for f_spm and f_offset

        elif filtyp == 'FRFP':
            outfile = out2 + '.P'
            print('\tSaving frequency residual fit parameters to:\n\t\t'
                    + outfile)
            file_object = open(outfile, 'wb')
            pickle.dump(contents, file_object)
            file_object.close()

        elif filtyp == 'PNFP':
            outfile = out2 + '.P'
            print('\tSaving power normalization fit parameters to:\n\t\t'
                    + outfile)
            file_object = open(outfile, 'wb')
            pickle.dump(contents, file_object)
            file_object.close()




    return None
