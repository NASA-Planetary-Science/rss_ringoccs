# get_CORSS_8001.py
# install the CORSS_8001 Cassini ring occultation RSS archive set from NASA PDS
# Destination directory: rss_ringoccs/data/CORSS_8001/
# Directory structure:
# CORSS_8001/
# ├── browse
# ├── catalog
# └── data
#     ├── Rev007
#     │   ├── Rev007E
#     │   │   ├── Rev007E_RSS_2005_123_K34_E
#     │   │   ├── Rev007E_RSS_2005_123_S43_E
#     │   │   ├── Rev007E_RSS_2005_123_X34_E
#     │   │   └── Rev007E_RSS_2005_123_X43_E


import os

data_dir = os.getcwd()+'/../data/'
dest_dir = data_dir+'CORSS_8001'

remote_file = 'https://pds-rings.seti.org/holdings/archives-volumes/CORSS_8xxx/CORSS_8001.tar.gz'
cmd = 'curl -fsSL '+remote_file +' | tar xzf  - -C '+data_dir

# print(dest_dir)
# print(cmd)

if os.path.isdir(dest_dir):
    print(dest_dir,'already exists')
else:
    print('executing\n',cmd)
    response = os.system(cmd)
    print('All done!')
