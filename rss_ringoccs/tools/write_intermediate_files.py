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
    2018 Sep 21 - jfong - update seqnum for new dates
    2018 Sep 26 - jfong - update print statement
"""

import os
from time import strftime
import pdb
import pickle
from .date_to_rev import date_to_rev
import numpy as np

def write_intermediate_files(year, doy, band, dsn, profdir, filtyp,
        contents):
    """
    Write intermediate files (frequency offset file, residual fit parameters
        file, power normalization fit parameters file) to directory
        ../output/RevXXX/P/RevXXXP_RSS_YEAR_DOY_BDD_P/
                RSS_YEAR_DOY_BDD_filtype_DATE_#.ext
        where XXX is the 3-digit rev number, P is the profile direction
        (I for ingress, E for egress), YEAR is the 4-digit year,
        DOY is the 3-digit day-of-year, B is the 1-letter band, DD is the
        2-digit dsn station number, filtyp is the abbreviation for the output
        file type, DATE is the current date in YYYYMMDD format, # is the
        number of file types written so far in current date, ext is the
        extension associated with filtyp.

    Args:
        year (str):
            year of event
        doy (str):
            day-of-year of event
        band (str):
            Name of the wavelength of transmission (S, X, or K)
        dsn (str):
            Deep Space Network ID of the station that observed the spacecraft
            for this event
        profdir (str):
            Profile direction of event.
        filtyp (str):
            Type of intermediate file to search for. 'FOF' for frequency
            offset file, 'FRFP' for frequency residual fit parameters,
            and 'PNFP' for power normalization fit parameters.
        contents (dict):
            Dictionary containing contents to be written to intermediate file.
    """
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
    
        dirstr = ('../output/Rev' + rev + '/' + dd + '/' + 'Rev' + rev + 
                pd2 + dd + '_' + filestr + '/')

        # Create output file name without file extension
        curday = strftime('%Y%m%d')

        out1 = dirstr + filestr + '_' + filtyp.upper() + '_' + curday

        # Check if directory exists, if not, create it
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
            # Create new directory and write first file of filtyp
            print('\tCreating directory:\n\t\t' + dirstr)
            os.system('[ ! -d ' + dirstr + ' ] && mkdir -p ' + dirstr)

            seq_num = '0001'
            out2 = out1 + '_' + seq_num

        if filtyp == 'FOF':
            outfile = out2 + '.TXT'
            f_spm = contents['f_spm']
            f_offset = contents['f_offset']
            print('\tSaving frequency offset calculations to:\n\t\t'
                    + '/'.join(outfile.split('/')[0:5]) + '/\n\t\t\t'
                        + outfile.split('/')[-1])
            np.savetxt(outfile, np.c_[f_spm, f_offset], fmt='%32.16F'*2)

        elif filtyp == 'FRFP':
            outfile = out2 + '.P'
            print('\tSaving frequency residual fit parameters to:\n\t\t'
                    + '/'.join(outfile.split('/')[0:5]) + '/\n\t\t\t'
                        + outfile.split('/')[-1])
            file_object = open(outfile, 'wb')
            pickle.dump(contents, file_object)
            file_object.close()

        elif filtyp == 'PNFP':
            outfile = out2 + '.P'
            print('\tSaving power normalization fit parameters to:\n\t\t'
                    + '/'.join(outfile.split('/')[0:5]) + '/\n\t\t\t'
                        + outfile.split('/')[-1])
            file_object = open(outfile, 'wb')
            pickle.dump(contents, file_object)
            file_object.close()

    return None
