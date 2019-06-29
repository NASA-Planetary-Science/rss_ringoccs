"""
Purpose: Read Uranus data files.

NOTE: Uranus data files are on our local TC2017/data directory!!
NOTE: 1 second gaps between egress UC0482 and UC0483 tapes
"""
import pdb
import numpy as np
from scipy.fftpack import hilbert
import os
import sys

class VGRUranusReader(object):
    def __init__(self,data_dir,prof_dir):
        # check if data directory exists
        print('Reading data files...')
        if not os.path.isdir(data_dir) :
            print('Directory '+data_dir+' not found.')
            sys.exit()
        # check if provided direction is allowed
        if prof_dir.upper() == 'I' :
            file_list = [data_dir + x for x in [
                    'SCRA_1066/UC0461/F00001.DAT',
                    'SCRA_1066/UC0462/F00001.DAT',
                    'SCRA_1067/UC0463/F00001.DAT']]
        elif prof_dir.upper() == 'E' :
            file_list = [data_dir + x for x in [
                    'SCRA_1070/UC0482/F00001.DAT',
                    'SCRA_1070/UC0483/F00001.DAT',
                    'SCRA_1070/UC0484/F00001.DAT',
                    'SCRA_1070/UC0485/F00001.DAT']]
        else:
            print('Profile direction \''+dir+'\' not recognized.')
            sys.exit()
        # check if files are in the appropriate directory
        file_check = 0
        for file in file_list:
            if not os.path.exists(file):
                print('File '+file+'\n was not found. Check input'
                      +'\ndata directory to make sure file is '
                      +'in the appropriate location')
                file += 1
        if file_check > 0:
            sys.exit()

        hdr_bytes       = 32
        nrecords        = 24000
        record_bytes    = 4390
        rec_hdr_bytes   = 80
        rec_data_bytes  = 4000
        rec_PPM_bytes   = 300
        rec_offset_bytes= 10
        nchan           = 4

        n_words_of_record_info = 16
        n_words_of_hdr = 40
        n_words_of_data = 2000
        n_words_of_ppm = 150
        n_words_of_operator_entry = 5

        n_bytes_of_record_info = n_words_of_record_info*2
        n_bytes_of_hdr = n_words_of_hdr*2
        rec_hdr = int(rec_hdr_bytes/2)
        n_bytes_of_data = n_words_of_data*2
        n_bytes_of_ppm = n_words_of_ppm*2
        n_bytes_of_operator_entry = n_words_of_operator_entry*2

        channel         = 2 # only A/D that has good data

        # read in binary data
        data = []
        for input_file in file_list:
            data_read = self.read_vgr_uranus(input_file)
            data.append(data_read)
        data_out = [item for sublist in data for item in sublist]

        # get I and Q
        I = np.array([x-128 for x in data_out])
        Q = np.array(hilbert(I))

        # generate SPM array
        Npts_per_file   = nrecords * rec_data_bytes/nchan
        nfiles = len(file_list)
        samples_per_sec = 50000.
        dt              = 1./samples_per_sec

        if prof_dir.upper() == 'I' :
            # 1986/024T22:27:01.000000
            spm0 = 22.*3600. + 27. * 60. + 1.0
            spm  = spm0 + dt*np.arange(Npts_per_file*nfiles)
        elif prof_dir.upper() == 'E' :
            # 'SPM (1986 Jan 25)' ; notice 1 sec gap between tapes
            spm0 =  1.*3600. + 15. * 60. + 2. # 1986/025T01:15:02.000000
            spm1 =  1.*3600. + 23. * 60. + 3. # 1986/025T01:23:03.000000
            spm2 =  1.*3600. + 31. * 60. + 3. # 1986/025T01:31:03.000000
            spm3 =  1.*3600. + 39. * 60. + 3. # 1986/025T01:39:03.000000
            spm  = np.reshape([[spm0 + dt*np.arange(Npts_per_file)]
                            + [spm1 + dt*np.arange(Npts_per_file)]
                            + [spm2 + dt*np.arange(Npts_per_file)]
                            + [spm3 + dt*np.arange(Npts_per_file)]],
                            len(I))

        self.spm = spm
        self.I = I
        self.Q = Q
        self.prof_dir = prof_dir


    def read_vgr_uranus(self,input_file):

        hdr_bytes       = 32
        nrecords        = 24000
        record_bytes    = 4390
        rec_hdr_bytes   = 80
        rec_data_bytes  = 4000
        rec_PPM_bytes   = 300
        rec_offset_bytes= 10
        nchan           = 4

        n_words_of_record_info = 16
        n_words_of_hdr = 40
        n_words_of_data = 2000
        n_words_of_ppm = 150
        n_words_of_operator_entry = 5

        n_bytes_of_record_info = n_words_of_record_info*2
        n_bytes_of_hdr = n_words_of_hdr*2
        rec_hdr = int(rec_hdr_bytes/2)
        n_bytes_of_data = n_words_of_data*2
        n_bytes_of_ppm = n_words_of_ppm*2
        n_bytes_of_operator_entry = n_words_of_operator_entry*2

        channel         = 2 # only A/D that has good data

        with open(input_file, 'rb') as f:
            rec_info = f.read(n_bytes_of_record_info)
            data_raw_full = []
            for n in range(nrecords):
                hdr_raw = f.read(n_bytes_of_hdr)
                data_raw = f.read(n_bytes_of_data)
                data_raw_full.append(data_raw)
                ppm_raw = f.read(n_bytes_of_ppm)
                operator_entry_raw = f.read(n_bytes_of_operator_entry)
            data_out = [item for sublist in data_raw_full for item in sublist]
            data = np.reshape(data_out, (int(nrecords*rec_data_bytes/nchan),nchan))[:,channel]
            return data

