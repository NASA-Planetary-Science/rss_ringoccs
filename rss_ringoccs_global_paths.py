# rss_ringoccs_global_paths.py

# Edit the following line with the full absolute path to your rss_ringoccs/ directory
# Be sure to include the ending slash.

global_path_to_rss_ringoccs='/Volumes/dione_raid2/Research/PDART_RSS/rss_ringoccs/'

# edit the following line with  the full absolute path to your rss_ringoccs_local/ directory
# that contains routines that call rss_ringoccs.
# Be sure to include the ending slash

global_path_to_rss_ringoccs_local='/Volumes/dione_raid2/Research/PDART_RSS/rss_ringoccs_local/'


'''
Usage:
To import these definitions into program.py located in otherdir/:

cd otherdir/
ln -s <path_to_rss_ringoccs>/rss_ringoccs_global_paths.py

in program.py:

import rss_ringoccs_global_paths # required to allow program to run from other directories

global global_path_to_rss_ringoccs
global_path_to_rss_ringoccs = rss_ringoccs_global_paths.global_path_to_rss_ringoccs

global global_path_to_local # note shortened name in next line
global_path_to_local = rss_ringoccs_global_paths.global_path_to_rss_ringoccs_local

 Add any desired auxiliary path definitions based on above two globals
        
global global_path_to_data 
global_path_to_data = global_path_to_rss_ringoccs + 'data/' 
        
global path_to_local_output
global_path_to_local_output = global_path_to_local + 'output/'
...

Then, for example, within the program:

geofile = global_path_to_data + 'Rev007/Rev007E/Rev007E_RSS_2005_123_X43_E/RSS_2005_123_X43_E_GEO_20250408_0001.TAB'
myoutputfile = global_path_to_output+ 'myoutputfilename.txt'
'''

