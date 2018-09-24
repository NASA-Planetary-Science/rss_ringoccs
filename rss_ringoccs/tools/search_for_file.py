"""

search_for_file.py

Purpose: This function searches for an existing frequency offset text file,
         frequency offset fit parameter pickle file, or a power normalization
         fit parameter pickle file. If more than one of the desired files 
         exist, then the most recent file will be used.

Revisions:
    2018 Sep 17 - jfong - original
    2018 Sep 24 - jfong - hardcode relative path (../ouput/*)
"""

import glob
import os
import pdb
import sys
from .date_to_rev import date_to_rev

def search_for_file(year, doy, band, dsn, profdir, filtyp):
    # Remove 'DSS-' from dsn name
    dsn = dsn.split('-')[-1]

    # Rename profdir to be I for ingress, E for egress, CI/CE for chord
    pd1 = (profdir.split('"')[1])[0]
    if pd1 == 'B':
        pd1 = '*'

    # Retrieve rev number in 3-digit format
    rev = date_to_rev(year, doy) 
    
    # Construct directory structure nomenclature
    #   RevXXX/P/RevXXXP_RSS_YEAR_DOY_BDD_P/RSS_YEAR_DOY_BDD_filtype_DATE_#
    # sample chord: Rev054/E/Rev054CE_RSS_2007_353_K55_E/
    # Rev054/E/Rev054CE_RSS_2007_353_K55_E
#    pdb.set_trace()
    filestr = ('RSS_' + str(year) + '_' + str(doy) + '_' + str(band) +
            str(dsn) + '_' + pd1)
    dirstr = ('../output/Rev' + rev + '/' + pd1 + '/' + 'Rev' + rev + pd1 +
            '_' + filestr + '/')
    
#    fsearch = dirstr + filestr + '*' + filtyp.upper() + '*' + '.*'
    #cwds = (os.getcwd()).split('/')
    # index of first appearance of rss_ringoccs
    #ind = cwds.index('rss_ringoccs')
    #relpath = '../' * int(len(cwds)-ind-1)
    
    #fsrch = relpath + dirstr +  '*' + filtyp.upper() + '*' + '.*'
    fsrch = dirstr +  '*' + filtyp.upper() + '*' + '.*'

    all_files = []
    for filename in glob.glob(fsrch, recursive=True):
        all_files.append(filename)
    

    if len(all_files) > 1:
        # sort by date, oldest to newest
        datestr = [x.split('/')[-1][-12:-4] for x in all_files]
        all_files_sorted = [x for y, x in sorted(zip(datestr,all_files))]
        rfile = all_files_sorted[-1]
    elif len(all_files) == 0:
        print('WARNING (search_for_file.py): '+ filtyp + ' not found!')
        rfile = 'N/A'
    else:
        rfile = all_files[0]


#    search_dir = 'Rev' + rev + '/' +
    return rfile
