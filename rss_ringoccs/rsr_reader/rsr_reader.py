#!/usr/bin/env python
"""

rsr_reader.py

Purpose: Class to create an instance linked to an RSR file.

Revisions:
      gjs_rsr_reader_v2.py
   2018 Feb 14 - gsteranka - Original version
   2018 Feb 20 - gsteranka - Made edits JWF suggested in her copy of the code.
                             Reading RSR header in __init__, since you always
                             have to do that anyway. Class name changed
                             from ReadRSR to RSRReader. "read_hdr" method now
                             private, since user doesn't need to use
   2018 Feb 20 - gsteranka - DSN output has "DSS-" in front of the number
      gjs_rsr_reader_v3.py
   2018 Feb 28 - gsteranka - Copied from v2. Edited names of methods and
                             edited docstrings. Also changed SPM_vals
                             definition so there's no rounding error
   2018 Mar 02 - gsteranka - Edited so predicted sky frequency and I/Q info
                             are not set as attributes to the objects, but
                             rather returned into variables. Changed so now you
                             can't make two different objects with different
                             attributes set. Also eliminated separate header
                             reading routine and am just doing it in __init__
   2018 Mar 07 - gsteranka - Fixed bug where length of SPM_vals outputted by
                             get_IQ doesn't match length of IQ_m
      rsr_reader.py
   2018 Mar 20 - gsteranka - Copy to official version and remove debug steps
   2018 Mar 21 - gsteranka - Added keyword argument in get_IQ method to
                             decimate 16kHz files to 1kHz sampling if True
   2018 May 08 - gsteranla - Added input checks
   2018 May 30 - gsteranka - Added history attribute

*************************VARIABLES*************************
NAME - TYPE - scalar/array - PURPOSE
endian - str - scalar - Endian format to use in reading rsr_file
f_spm - float - array - Frequency SPM to evaluate sky frequency over
f_sky_pred - float - array - Sky frequency read in from polynomials in RSR file
rsr_file - str - scalar - Full path name to an RSR file
sh_sample_rate_hz - int - scalar - Sample rate of RSR file in Hz
spm_range - float - array - Range of SPM to read data from RSR file
spm_vals - float - array - Seconds past midnight of data measurements
"""

import multiprocessing
from multiprocessing import Process
from multiprocessing import Queue
import numpy as np
import os
import platform
from scipy.signal import decimate
import struct
import sys
import time

class RSRReader(object):
    """
    Reads the header of a raw RSR file when you first create an instance.
    Then reads the full RSR file to calculate predicted sky frequency at the
    specified times, or reads the full RSR file to read the raw measured I
    and Q values.

    Example:
        >>> # Define instance and set header attributes
        >>> rsr_inst = RSRReader(rsr_file)
        >>>
        >>>  # Get predicted sky frequency at chosen SPM values f_spm
        >>> f_spm_returned, f_sky_pred = rsr_inst.get_f_sky_pred(f_spm=f_spm)
        >>>
        >>> # Get raw measured I and Q, and raw SPM, over specified spm_range
        >>> spm_vals, IQ_m = rsr_inst.get_IQ(spm_range=spm_range)

    Data attributes:
        rsr_file (str): Full path name of a raw RSR file to read
        spm_vals (np.ndarray): Seconds Past Midnight array of times over
            entire rsr file
        doy (int): Day of year of event
        year (int): Year of event
        dsn (str): Deep Space Network ID of the station that observed the
            spacecraft for this event
        band (bytes): Name of the wavelength of transmission (S, X, or Ka)
        sample_rate_khz (int): Sample rate, in kHz, of transmission (1 or 16)
        history (dict): Dictionary recording parameters of the run
    """

    # Define field names and format variables used in each method
    # (from here until __init__)
    __endian = '>'

    # Header contents of SFDU
    __sfdu_field_names = [
        'sfdu_authority1', 'sfdu_authority2', 'sfdu_authority3',
        'sfdu_authority4',
        'sfdu_version',
        'sfdu_class',
        'sfdu_reserved1', 'sfdu_reserved2',
        'sfdu_data_desription1', 'sfdu_data_desription2',
        'sfdu_data_desription3', 'sfdu_data_desription4',
        'sfdu_length']
    __sfdu_format = 'cccc'+'c'+'c'+'cc'+'cccc'+'Q'

    # Header Aggregation
    __ha_field_names = ['ha_type', 'ha_length']
    __ha_format = 'H'+'H'

    # Primary Header
    __ph_field_names = [
        'ph_type',
        'ph_length',
        'ph_data_major',
        'ph_data_minor',
        'ph_mission_ID',
        'ph_format_code']
    __ph_format = 'H'+'H'+'B'+'B'+'B'+'B'

    # Secondary Header
    __sh_field_names = [
        'sh_type',
        'sh_length',
        'sh_originator',
        'sh_last_modifier',
        'sh_rsr_software_id',
        'sh_record_sequence_number',
        'sh_spc_id',
        'sh_dss_id',
        'sh_rsr_id',
        'sh_schan_id',
        'sh_reserved',
        'sh_spacecraft',
        'sh_prdx_pass_number',
        'sh_ulband',
        'sh_dl_band',
        'sh_trk_mode',
        'sh_ul_dss_id',
        'sh_fgain_px_no',
        'sh_fgain_if_bandwidth',
        'sh_frov_flag',
        'sh_attenuation',
        'sh_adc_rms','sh_adc_peak',
        'sh_year','sh_doy','sh_seconds',
        'sh_bits_per_sample',
        'sh_data_error',
        'sh_sample_rate',
        'sh_ddc_lo',
        'sh_rfif_lo','sh_sfdu_year','sh_sfdu_doy','sh_sfdu_seconds',
        'sh_predicts_time_shift','sh_predicts_freq_override',
        'sh_predicts_freq_rate','sh_predicts_freq_offset',
        'sh_sub_channel_freq',
        'sh_rf_freq_point_1','sh_rf_freq_point_2','sh_rf_freq_point_3',
        'sh_schan_freq_point_1','sh_schan_freq_point_2',
        'sh_schan_freq_point_3',
        'sh_schan_freq_poly_coef_1','sh_schan_freq_poly_coef_2',
        'sh_schan_freq_poly_coef_3',
        'sh_schan_accum_phase',
        'sh_schan_phase_poly_coef_1','sh_schan_phase_poly_coef_2',
        'sh_schan_phase_poly_coef_3','sh_schan_phase_poly_coef_4',
        'sh_reserved2a','sh_reserved2b']
    __sh_format = 'hh'+'BBh'+'hBB'+'BBcBHccBBbBBBBBHHIBBHHHHH'+22*'d'

    # Data
    __data_field_names = [
        'Data_type',
        'Data_length',
        'Data_QI']
    __data_header_format = 'HH'

    __field_names = (__sfdu_field_names + __ha_field_names + __ph_field_names
        + __sh_field_names + __data_field_names)


    def __init__(self,rsr_file, decimate_16khz_to_1khz=False,
            cpu_count=multiprocessing.cpu_count(), TEST=False):
        """
        Purpose:
        Sets full path name of RSR file as an attribute to the instance, and
        reads the header of the RSR file, which sets the RSR header attributes
        spm_vals, doy, year, dsn, band, sample_rate_khz. Also sets private
        attributes for reading the RSR file later in the "get" methods.

        Args:
            rsr_file (str): Full path name of a raw RSR file to read
            decimate_16khz_to_1khz (bool): Optional Boolean argument which, if
                set to True, decimates 16kHz files down to 1kHz sampling rate.
                Note that this is a sticky keyword - if you set it to True, it
                will be True for any subsequent calls from the instance until
                you explicitly set it to False. This keyword is linked to the
                private attribute __decimate_16khz_to_1khz
            cpu_count (int): Number of cores to use when reading data in from
                file. Default is number of cores on your computer
            TEST (bool): Optional boolean variable which, when set to True,
                prints the header attributes that were set

        Dependencies:
            [1] RSRReader
            [2] numpy
            [3] os
            [4] platform
            [5] scipy.signal.decimate
            [6] struct
            [7] sys
            [8] time

        Warnings:
            [1] If you try setting decimate_16khz_to_1khz=True for a 1kHz file,
                it will just ignore you
            [2] 16kHz files will take a few minutes to read
        """

        # Ensure TEST is Boolean
        if type(TEST) != bool:
            print('WARNING (RSRReader): TEST input should be Boolean. Assuming '+
                'False')
            TEST = False

        if type(cpu_count) != int:
            print('WARNING (RSRReader): cpu_count keyword should be an '
                + 'integer. Setting to number of cores on your computer.')
            cpu_count = multiprocessing.cpu_count()

        self.rsr_file = rsr_file

        # Default argment for __set_IQ and cpu_count
        self.__decimate_16khz_to_1khz = decimate_16khz_to_1khz
        self.__cpu_count = cpu_count

        # Record information about the run
        self.__set_history()

        # Length of SFDU header, and a function to read the header in the
        # proper format
        struct_hdr_fmt = (self.__endian + self.__sfdu_format + self.__ha_format
            + self.__ph_format + self.__sh_format
            + self.__data_header_format)
        struct_hdr_len = struct.calcsize(struct_hdr_fmt)
        struct_unpack_hdr = struct.Struct(struct_hdr_fmt).unpack_from

        # Open first header of SFDU and put into format for unpacking below
        try:
            with open(self.rsr_file,'rb') as f:
                sfdu_hdr_raw = f.read(struct_hdr_len)
                f.close()
        except FileNotFoundError as err:
            print('ERROR (RSRReader): File not found! {}'.format(err))
            sys.exit()

        # Unpack SFDU header
        sfdu_hdr = struct_unpack_hdr(sfdu_hdr_raw)

        # Able to use field names to reference parts of SFDU header
        sfdu_hdr_dict = dict(zip(self.__field_names,sfdu_hdr))

        # Find number of SFDU in file, and number of points per SFDU
        rsr_size = os.path.getsize(self.rsr_file)
        bytes_per_sfdu = sfdu_hdr_dict['sfdu_length'] + 20
        n_sfdu = rsr_size / bytes_per_sfdu
        sh_bits_per_sample = sfdu_hdr_dict['sh_bits_per_sample']
        bytes_per_sample = sh_bits_per_sample / 8
        data_length_per_sfdu = sfdu_hdr_dict['Data_length']
        n_pts_per_sfdu = np.int(data_length_per_sfdu / (2*bytes_per_sample))

        # Get array of SPM values for whole file
        sh_sample_rate_hz = sfdu_hdr_dict['sh_sample_rate'] * 1000.0
        sh_sfdu_seconds = sfdu_hdr_dict['sh_sfdu_seconds']
        dt = 1.0 / sh_sample_rate_hz
        end_spm_of_rsr = sh_sfdu_seconds + n_pts_per_sfdu*n_sfdu*dt
        n_pts = round((end_spm_of_rsr - sh_sfdu_seconds)/dt)
        spm_vals = float(sh_sfdu_seconds) + dt*np.arange(n_pts)

        # Set RSR header attributes
        self.spm_vals = spm_vals
        self.doy = sfdu_hdr_dict['sh_doy']
        self.year = sfdu_hdr_dict['sh_year']
        self.dsn = 'DSS-'+str(sfdu_hdr_dict['sh_dss_id'])
        self.band = sfdu_hdr_dict['sh_dl_band']
        self.sample_rate_khz = sfdu_hdr_dict['sh_sample_rate']

        # Set attributes for later reading of rest of RSR file
        self.__n_pts_per_sfdu = n_pts_per_sfdu
        self.__n_sfdu = n_sfdu

        if TEST:
            print('First 10 raw SPM:')
            print(self.spm_vals[0:10])
            print('Year, DOY, DSN, band, sample_rate:')
            print(str(self.year) + ', ' + str(self.doy) + ', '
                + str(self.dsn) + ', ' + str(self.band) + ', '
                + str(self.sample_rate_khz))

        self.__set_IQ(TEST=TEST)


    def __set_sfdu_unpack(self, spm_range):
        """
        Set private attribute __sfdu_unpack, which is used to unpack the
        RSR file one SFDU at a time. Also sets attributes for the start and
        end SFDU to read. Not included in __init__ because it's not necessary
        for reading the header information
        
        Args:
            spm_range (list): 2-element array of range of SPM values to read
                over. Passed from either set_f_sky_pred or set_IQ
        """

        # Specify which SFDUs you want to read
        start_spm_ind = np.argmin(abs(self.spm_vals - spm_range[0]))
        end_spm_ind = np.argmin(abs(self.spm_vals - spm_range[1]))
        start_sfdu = int(start_spm_ind / self.__n_pts_per_sfdu)
        end_sfdu = int(end_spm_ind / self.__n_pts_per_sfdu)
        if end_sfdu > self.__n_sfdu:
            end_sfdu = self.__n_sfdu
        spm_vals = self.spm_vals[self.__n_pts_per_sfdu*start_sfdu:
            self.__n_pts_per_sfdu*end_sfdu]

        # Format in which to read rest of RSR file one SFDU at a time
        data_format = (self.__data_header_format
            + np.int(self.__n_pts_per_sfdu)*'hh')

        # Format to read RSR file in
        rsr_fmt = (self.__endian + self.__sfdu_format + self.__ha_format
            + self.__ph_format + self.__sh_format + data_format)
        sfdu_len = struct.calcsize(rsr_fmt)
        sfdu_unpack = struct.Struct(rsr_fmt).unpack_from

        # Define structure of RSR
        with open(self.rsr_file,'rb') as f:
            rsr_struct = f.read()
            f.close()

        # Define private attributes of object
        self.__sfdu_unpack = sfdu_unpack
        self.__rsr_struct = rsr_struct
        self.__start_sfdu = start_sfdu
        self.__end_sfdu = end_sfdu
        self.__sfdu_len = sfdu_len


    def get_f_sky_pred(self, f_spm=None, TEST=False):
        """
        Calculate predicted sky frequency at user-defined times using
        polynomial coefficients in each SFDU. Returns f_spm and
        f_sky_pred
        
        Args:
            f_spm (np.ndarray): Array of SPM values to evaluate predicted
                sky frequency at. Default is at 1 second spacing over entire
                data set.
            TEST (bool): Print the first few predicted sky frequency values if
                set to True

        Outputs:
            f_spm (np.ndarray): Array of SPM values that predicted sky frequency
                was evaluated at.
            f_sky_pred (np.ndarray): Predicted sky frequency, calculated from
                the polynomial coefficients in the RSR file

        Dependencies:
            [1] numpy
            [2] struct
            [3] sys

        Warnings:
            [1] Will take a few minutes to run for 16kHz files
        """

        # Ensure TEST is Boolean
        if type(TEST) != bool:
            print('WARNING (RSRReader): TEST input should be Boolean. Assuming '
                + 'False')
            TEST = False

        spm_vals = self.spm_vals

        # Default 1 second spacing over range of spm_vals
        if f_spm is None:
            f_spm = np.arange(min(spm_vals),max(spm_vals),1.0)

        # Input f_spm needs to be a numpy array
        if type(f_spm) != np.ndarray:
            print('ERROR (RSRReader.get_f_sky_pred): f_spm must be a numpy '
                + 'array')
            sys.exit()

        # SPM is often just a little off from what it's supposed to be, so need
        #     to later round to a certain number of decimals. 16 kHz files have
        #     more decimals
        if self.sample_rate_khz == 1:
            round_decimal = 4
        else:
            round_decimal = 7

        # Set up unpacking of full RSR. Catch error if f_spm is not array of
        #     floats/integers
        try:
            spm_range = [min(f_spm),max(f_spm)]
            self.__set_sfdu_unpack(spm_range)
        except TypeError as err:
            print('f_spm must be an array of floats or integers: '+
                '{}'.format(err))
            sys.exit()

        # Arrays to contain sets of frequency polynomials
        rfif_lo_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        ddc_lo_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        freq_poly1_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        freq_poly2_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        freq_poly3_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        time_stamp_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        n_iter = 0
        for i_sfdu in range(self.__start_sfdu,self.__end_sfdu+1):
            sfdu = self.__rsr_struct[i_sfdu*self.__sfdu_len:
                i_sfdu*self.__sfdu_len + self.__sfdu_len]

            # If end of file reached
            if len(sfdu) == 0:
                break

            # Unpack SFDU into readable format
            s = self.__sfdu_unpack(sfdu)
            s_dict = dict(zip(self.__field_names, s))

            # Beginning of time range for frequency polynomials for this SFDU
            _time_stamp = spm_vals[i_sfdu*self.__n_pts_per_sfdu]

            rfif_lo_array[n_iter] = s_dict['sh_rfif_lo']
            ddc_lo_array[n_iter] = s_dict['sh_ddc_lo']
            freq_poly1_array[n_iter] = s_dict['sh_schan_freq_poly_coef_1']
            freq_poly2_array[n_iter] = s_dict['sh_schan_freq_poly_coef_2']
            freq_poly3_array[n_iter] = s_dict['sh_schan_freq_poly_coef_3']
            time_stamp_array[n_iter] = round(_time_stamp,round_decimal)

            n_iter += 1

        # Calculate sky frequency for the input time array
        f_sky_pred = np.zeros(len(f_spm))
        for i in range(len(f_spm)):
            # Find correct frequency info for time
            _ind = np.arange(len(time_stamp_array))
            _time_ind = (_ind[f_spm[i] >= time_stamp_array])[-1]

            _msec = f_spm[i] - f_spm[i].astype(int)
            f_sky_pred[i] = ((rfif_lo_array[_time_ind]
                + ddc_lo_array[_time_ind])*1.0e6
                - freq_poly1_array[_time_ind]
                - freq_poly2_array[_time_ind]*_msec
                - freq_poly3_array[_time_ind]*(_msec**2))

        if TEST:
            for i in range(10):
                print('%24.16f    %30.16f' % (f_spm[i],f_sky_pred[i]))

        # Return f_spm and evaluated predicted sky frequency
        return f_spm, f_sky_pred


    def __set_IQ(self, TEST=False):
        """
        Read full RSR file to find the raw measured I and Q over the
        specified spm_range of the file. Returns raw SPM values over the
        specified time range and the raw measured I and Q values.

        Args:
            TEST (bool): If True, print first 10 raw measured I and Q

        Adds attributes:
            spm_vals (np.ndarray): Raw resolution SPM values over specified
                spm_range
            IQ_m (np.ndarray): Raw measured complex signal over the specified
                spm_range
        """

        # Ensure TEST is Boolean
        if type(TEST) != bool:
            print('WARNING (RSRReader): TEST input should be Boolean. Assuming '
                + 'False')
            TEST = False

        spm_vals = self.spm_vals

        # Ensure that input is a Boolean
        if type(self.__decimate_16khz_to_1khz) != bool:
            print('WARNING (RSRReader.get_IQ): Expected Boolean input for '+
                'decimate_16khz_to_1khz keyword. Ignoring input')
            self.__decimate_16khz_to_1khz = False

        decimate_16khz_to_1khz = self.__decimate_16khz_to_1khz

        spm_range = [min(spm_vals), max(spm_vals)]
        self.__set_sfdu_unpack(spm_range)

        # Reduce SPM array to match the I and Q arrays to be made
        spm_vals = self.spm_vals[self.__n_pts_per_sfdu*self.__start_sfdu:
            self.__n_pts_per_sfdu*(self.__end_sfdu+1)]

        # Multiprocessing to retrieve data from RSR file
        results = []
        queues = [Queue() for i in range(self.__cpu_count)]
        n_loops = self.__end_sfdu - self.__start_sfdu + 1
        n_per_core = int(np.floor(n_loops/self.__cpu_count))
        loop_args = [(i*n_per_core, (i+1)*n_per_core, n_loops,
            queues[i]) for i in range(self.__cpu_count)]
        loop_args[-1] = ((self.__cpu_count-1)*n_per_core,
            self.__end_sfdu + 1, n_loops, queues[-1])
        jobs = [Process(target=self.__loop, args=(a)) for a in loop_args]
        for j in jobs: j.start()
        for q in queues: results.append(q.get())
        for j in jobs: j.join()
        IQ_m = np.hstack(results)
        print('\n')

        # Decimate 16kHz file to 1kHz spacing if specified
        if decimate_16khz_to_1khz & (self.sample_rate_khz == 16):
            IQ_m = decimate(IQ_m, 4, zero_phase=True)
            IQ_m = decimate(IQ_m, 4, zero_phase=True)

            n_pts = len(IQ_m)
            dt = 1.0/float(1000)
            spm_vals = spm_vals[0] + dt*np.arange(n_pts)
        elif decimate_16khz_to_1khz & (self.sample_rate_khz == 1):
            print('WARNING (RSRReader.get_IQ): Cannot decimate a 1 kHz file '
                + 'any further. Skipping extra decimation')

        if TEST:
            print('First 10 SPM, I, and Q:')
            for i in range(10):
                print('%24.16f %15i %15i' %
                    (spm_vals[i], I_m[i], Q_m[i]))

        # Set spm_vals attribute corresponding to speified spm_range and the raw
        # measured I and Q
        self.spm_vals = spm_vals
        self.IQ_m = IQ_m


    def __loop(self, i_start, i_end, n_loops, queue=0):
        """
        Function to perform loop for multiprocessing
        """

        I_array = np.zeros(len(range(i_start, i_end))*self.__n_pts_per_sfdu)
        Q_array = np.zeros(len(range(i_start, i_end))*self.__n_pts_per_sfdu)
        i_iter = 0
        for i_sfdu in range(i_start, i_end):
            sfdu = self.__rsr_struct[i_sfdu*self.__sfdu_len:
                i_sfdu*self.__sfdu_len + self.__sfdu_len]

            # If EOF is reached
            if len(sfdu) == 0:
                break

            # Unpack SFDU into readable format
            s = self.__sfdu_unpack(sfdu)
            s_dict = dict(zip(self.__field_names, s))

            I_array[i_iter*self.__n_pts_per_sfdu:
                (i_iter+1)*self.__n_pts_per_sfdu] = (
                s[-2*self.__n_pts_per_sfdu+1::2])
            Q_array[i_iter*self.__n_pts_per_sfdu:
                (i_iter+1)*self.__n_pts_per_sfdu] = (
                s[-2*self.__n_pts_per_sfdu::2])

            i_iter += 1

        queue.put(I_array + 1j*Q_array)


    def __set_history(self):
        """
        Set Python dictionary recording information about the run as the
        history attribute
        """
        input_var_dict = {'rsr_file': self.rsr_file}
        input_kw_dict = {
            'decimate_16khz_to_1khz': self.__decimate_16khz_to_1khz}
        history_dict = {'User Name': os.getlogin(),
            'Host Name': os.uname().nodename,
            'Run Date': time.ctime() + ' ' + time.tzname[0],
            'Python Version': platform.python_version(),
            'Operating System': os.uname().sysname,
            'Source File': __file__.split('/')[-1],
            'Source Directory': __file__.rsplit('/',1)[0] +'/',
            'Input Variables': input_var_dict,
            'Input Keywords':input_kw_dict}
        self.history = history_dict


if __name__ == '__main__':

    rsr_file = '../../../../../data/s10-rev07-rsr-data/S10EAOE2005_123_0740NNNX43D.2A1'
    #rsr_file = '../../../../../data/s10-rev07-rsr-data/s10sroe2005123_0740nnnx43rd.2a2'
    rsr_inst = RSRReader(rsr_file)
