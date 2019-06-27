#!/usr/bin/env python
"""

:Purpose:
    Class to create an instance linked to an RSR file

:Dependencies:
    #. multiprocessing
    #. numpy
    #. os
    #. scipy
    #. struct
    #. sys
    #. time
"""

import multiprocessing
from multiprocessing import Process
from multiprocessing import Queue
import numpy as np
import os
import pdb
from scipy.signal import decimate
import struct
import sys
import time

from ..tools.history import get_rev_info
from ..tools.history import write_history_dict


class RSRReader(object):
    """
    :Purpose:
    Reads the header of a raw RSR file when you first create an instance.
    Then reads the full RSR file to read in the raw measured complex signal :math:`I+iQ`

    Arguments
        :rsr_file (*str*):
            Full path name of a raw RSR file to read. RSR files
            can be downloaded using the shell script in the data"
            directory of the GitHub clone

    Keyword Arguments
        :decimate_16khz_to_1khz (*bool*):
            Optional Boolean argument which, if
            set to True, decimates 16kHz files down to 1kHz sampling rate.
            Note that this is a sticky keyword - if you set it to True, it
            will be True for any subsequent calls from the instance until
            you explicitly set it to False. This keyword is linked to the
            private attribute __decimate_16khz_to_1khz
        :cpu_count (*int*):
            Number of cores to use when reading data in from
            file. Default is number of cores on your computer
        :verbose (*bool*):
            Optional boolean variable which, when set to True,
            prints the header attributes that were set

    Attributes
        :rsr_file (*str*): Full path name of a raw RSR file to read
        :spm_vals (*np.ndarray*): Seconds Past Midnight array of times over
            entire rsr file
        :doy (*int*): Day of year of event
        :year (*int*): Year of event
        :dsn (*str*): Deep Space Network ID of the downlink station
        :band (*str*): Name of the wavelength of downlink transmission
                        (S, X, or Ka)
        :ul_band (*str*): Name of the wavelength of uplink transmission
        :ul_dsn (*str*): Deep Space Network ID of uplink station
        :sample_rate_khz (*int*): Sample rate, in kHz, of transmission (1 or 16)
        :history (*dict*): Dictionary recording parameters of the run

    Example
        >>> # Import rss_ringoccs
        >>> import rss_ringoccs as rss
        >>> # Define instance and set header attributes, and read in raw data
        >>> rsr_inst = rss.rsr_reader.RSRReader(rsr_file)
        >>>  # Get predicted sky frequency at chosen SPM values f_spm
        >>> f_spm_returned, f_sky_pred = rsr_inst.get_f_sky_pred(f_spm=f_spm)

    Notes:
        #. Setting ``decimate_16khz_to_1khz=True`` for a 1kHz file will be ignored
        #. 16kHz files will take a few minutes to read and decimate
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
    __sfdu_format = 'cccc' + 'c' + 'c' + 'cc' + 'cccc' + 'Q'

    # Header Aggregation
    __ha_field_names = ['ha_type', 'ha_length']
    __ha_format = 'H' + 'H'

    # Primary Header
    __ph_field_names = [
        'ph_type',
        'ph_length',
        'ph_data_major',
        'ph_data_minor',
        'ph_mission_ID',
        'ph_format_code']
    __ph_format = 'H' + 'H' + 'B' + 'B' + 'B' + 'B'

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
        'sh_adc_rms', 'sh_adc_peak',
        'sh_year', 'sh_doy', 'sh_seconds',
        'sh_bits_per_sample',
        'sh_data_error',
        'sh_sample_rate',
        'sh_ddc_lo',
        'sh_rfif_lo', 'sh_sfdu_year', 'sh_sfdu_doy', 'sh_sfdu_seconds',
        'sh_predicts_time_shift', 'sh_predicts_freq_override',
        'sh_predicts_freq_rate', 'sh_predicts_freq_offset',
        'sh_sub_channel_freq',
        'sh_rf_freq_point_1', 'sh_rf_freq_point_2', 'sh_rf_freq_point_3',
        'sh_schan_freq_point_1', 'sh_schan_freq_point_2',
        'sh_schan_freq_point_3',
        'sh_schan_freq_poly_coef_1', 'sh_schan_freq_poly_coef_2',
        'sh_schan_freq_poly_coef_3',
        'sh_schan_accum_phase',
        'sh_schan_phase_poly_coef_1', 'sh_schan_phase_poly_coef_2',
        'sh_schan_phase_poly_coef_3', 'sh_schan_phase_poly_coef_4',
        'sh_reserved2a', 'sh_reserved2b']
    __sh_format = 'hh' + 'BBh' + 'hBB' + 'BBcBHccBBbBBBBBHHIBBHHHHH' + 22 * 'd'

    # Data
    __data_field_names = [
        'Data_type',
        'Data_length',
        'Data_QI']
    __data_header_format = 'HH'

    __field_names = (__sfdu_field_names + __ha_field_names + __ph_field_names
        + __sh_field_names + __data_field_names)

    def __init__(self, rsr_file, decimate_16khz_to_1khz=True, verbose=False):

        if not isinstance(verbose, bool):
            print('WARNING (RSRReader): verbose input should be Boolean. '
                + 'Assuming False. If you\'re trying to use 1 or 0, then you '
                +'should use the built-in Python booleans instead')
            verbose = False

        if not isinstance(decimate_16khz_to_1khz, bool):
            print('WARNING (RSRReader.get_IQ): Expected Boolean input for '
                + 'decimate_16khz_to_1khz keyword. Ignoring input. If you\'re '
                + 'trying to use 1 or 0, then you should use the built-in '
                + 'Python booleans instead')
            decimate_16khz_to_1khz = False

        self.rsr_file = rsr_file

        # Default argment for __set_IQ and cpu_count
        self.__decimate_16khz_to_1khz = decimate_16khz_to_1khz
        self.__cpu_count = multiprocessing.cpu_count()

        # Record information about the run
        self.__set_history()

        if verbose:
            print('\nExtracting information from RSR file...')
            print('\tReading header information...')
        # Length of SFDU header, and a function to read the header in the
        # proper format
        struct_hdr_fmt = (self.__endian + self.__sfdu_format + self.__ha_format
            + self.__ph_format + self.__sh_format
            + self.__data_header_format)
        struct_hdr_len = struct.calcsize(struct_hdr_fmt)
        struct_unpack_hdr = struct.Struct(struct_hdr_fmt).unpack_from

        # Open first header of SFDU and put into format for unpacking below
        try:
            with open(self.rsr_file, 'rb') as f:
                sfdu_hdr_raw = f.read(struct_hdr_len)
                f.close()
        except FileNotFoundError as err:
            print('ERROR (RSRReader): File not found! {}'.format(err))
            sys.exit()

        # Unpack SFDU header
        sfdu_hdr = struct_unpack_hdr(sfdu_hdr_raw)

        # Able to use field names to reference parts of SFDU header
        sfdu_hdr_dict = dict(zip(self.__field_names, sfdu_hdr))

        # Find number of SFDU in file, and number of points per SFDU
        rsr_size = os.path.getsize(self.rsr_file)
        bytes_per_sfdu = sfdu_hdr_dict['sfdu_length'] + 20
        n_sfdu = int(rsr_size / bytes_per_sfdu)
        if rsr_size % bytes_per_sfdu != 0 and verbose is True:
            print('WARNING (RSRReader): file size not the same as expected!\n'
                    + '\t Using n_sfdu=' + str(n_sfdu)
                    + ' instead of n_sfdu=' + str(rsr_size/bytes_per_sfdu))
        sh_bits_per_sample = sfdu_hdr_dict['sh_bits_per_sample']
        bytes_per_sample = sh_bits_per_sample / 8
        data_length_per_sfdu = sfdu_hdr_dict['Data_length']
        n_pts_per_sfdu = np.int(data_length_per_sfdu / (2 * bytes_per_sample))

        # Set fgain for threshold optical depth calculation
        self.fxgain_px_no = sfdu_hdr_dict['sh_fgain_px_no']
        self.fgain_if_bandwidth = sfdu_hdr_dict['sh_fgain_if_bandwidth']

        # Get array of SPM values for whole file
        sh_sample_rate_hz = sfdu_hdr_dict['sh_sample_rate'] * 1000.0
        sh_sfdu_seconds = sfdu_hdr_dict['sh_sfdu_seconds']
        dt = 1.0 / sh_sample_rate_hz
        end_spm_of_rsr = sh_sfdu_seconds + n_pts_per_sfdu * n_sfdu * dt
        n_pts = round((end_spm_of_rsr - sh_sfdu_seconds) / dt)
        spm_vals = float(sh_sfdu_seconds) + dt * np.arange(n_pts)

        # Set RSR header attributes
        self.spm_vals = spm_vals
        self.doy = sfdu_hdr_dict['sh_doy']
        self.year = sfdu_hdr_dict['sh_year']
        self.dsn = 'DSS-' + str(sfdu_hdr_dict['sh_dss_id'])
        self.ul_dsn = 'DSS-'+rsr_file.split('/')[-1].split('_')[1][5:7]
        if self.ul_dsn == 'DSS-MM':
            self.ul_dsn = 'DSS-'+str(sfdu_hdr_dict['sh_ul_dss_id'])
        self.band = chr(sfdu_hdr_dict['sh_dl_band'][0])
        self.ul_band = chr(sfdu_hdr_dict['sh_ulband'][0])
        # correct header info if not accurate
        if self.year > 2011 :
            self.track_mode = 2
            if int(sfdu_hdr_dict['sh_ul_dss_id']) < 10 :
                self.ul_dsn = 'DSS-'+rsr_file[-11:-9]
        else:
            self.track_mode = 1
        #print(self.ul_dsn)

        self.sample_rate_khz = sfdu_hdr_dict['sh_sample_rate']

        # DSS-74 files have the regular year and DOY set to 0
        if self.year == 0:
            self.year = sfdu_hdr_dict['sh_sfdu_year']
        if self.doy == 0:
            self.doy = sfdu_hdr_dict['sh_sfdu_doy']

        # Set attributes for later reading of rest of RSR file
        self.__n_pts_per_sfdu = n_pts_per_sfdu
        self.__n_sfdu = n_sfdu


        self.__set_IQ(verbose=verbose)

        # Set rev info for file creation
        self.__set_rev_info()
        if verbose:
            print('\t\tRev:\t\t\t' + self.rev_info['rev_num'])
            print('\t\tYear:\t\t\t' + str(self.year))
            print('\t\tDOY:\t\t\t' + str(self.doy))
            print('\t\tSPM range:\t\t' + str(self.spm_vals[0]) + ', '
                    + str(self.spm_vals[-1]))
            print('\t\tDSN:\t\t\t' + str(self.dsn))
            print('\t\tBand:\t\t\t' + str(self.band))
            print('\t\tSampling rate in kHz:\t' + str(self.sample_rate_khz))


    def __set_sfdu_unpack(self, spm_range):
        """
        Set private attribute ``__sfdu_unpack``, which is used to unpack the
        RSR file one SFDU at a time. Also sets attributes for the start and
        end SFDU to read. Not included in ``__init__`` because it's 
        not necessary for reading the header information

        Arguments
            :spm_range (*list*):
                2-element array of range of SPM values to read
                over. Passed from either ``set_f_sky_pred`` or ``set_IQ``
        """

        # Specify which SFDUs you want to read
        if self.sample_rate_khz == 1:
            start_spm_ind = np.argmin(abs(self.spm_vals - spm_range[0]))
            end_spm_ind = np.argmin(abs(self.spm_vals - spm_range[1]))
        elif self.sample_rate_khz == 16:
            start_spm_ind = np.argmin(abs(self.__spm_16khz - spm_range[0]))
            end_spm_ind = np.argmin(abs(self.__spm_16khz - spm_range[1]))
        start_sfdu = int(start_spm_ind / self.__n_pts_per_sfdu)
        end_sfdu = int(end_spm_ind / self.__n_pts_per_sfdu)
        if end_sfdu > self.__n_sfdu:
            end_sfdu = self.__n_sfdu

        # Format in which to read rest of RSR file one SFDU at a time
        data_format = (self.__data_header_format
            + np.int(self.__n_pts_per_sfdu) * 'hh')

        # Format to read RSR file in
        rsr_fmt = (self.__endian + self.__sfdu_format + self.__ha_format
            + self.__ph_format + self.__sh_format + data_format)
        sfdu_len = struct.calcsize(rsr_fmt)
        sfdu_unpack = struct.Struct(rsr_fmt).unpack_from

        # Define structure of RSR
        with open(self.rsr_file, 'rb') as f:
            rsr_struct = f.read()
            f.close()

        # Define private attributes of object
        self.__sfdu_unpack = sfdu_unpack
        self.__rsr_struct = rsr_struct
        self.__start_sfdu = start_sfdu
        self.__end_sfdu = end_sfdu
        self.__sfdu_len = sfdu_len

    def get_f_sky_pred(self, f_spm=None, verbose=False):
        """
        Calculate predicted sky frequency at user-defined times using
        polynomial coefficients in each SFDU.

        Arguments
            :f_spm (*np.ndarray*):
                Array of SPM values to evaluate predicted
                sky frequency at. Default is at 1 second spacing over entire
                data set.
            :verbose (*bool*):
                Print the first few predicted sky frequency values
                if set to True

        Returns
            :f_spm (*np.ndarray*):
                Array of SPM values that predicted sky
                frequency was evaluated at.
            :f_sky_pred (*np.ndarray*):
                Predicted sky frequency, calculated from
                the polynomial coefficients in the RSR file
        """

        if not isinstance(verbose, bool):
            print('WARNING (RSRReader): verbose input should be Boolean. '
                + 'Assuming False. If you\'re trying to use 1 or 0, then you '
                +'should use the built-in Python booleans instead')
            verbose = False

        if self.sample_rate_khz == 1:
            spm_vals = self.spm_vals
        elif self.sample_rate_khz == 16:
            spm_vals = self.__spm_16khz
        else:
            print('ERROR (RSRReader.f_sky_pred()): Sample rate must be either'
                + '16kHz or 1kHz')
            sys.exit()

        # Default 1 second spacing over range of spm_vals
        if f_spm is None:
            f_spm = np.arange(min(spm_vals), max(spm_vals), 1.0)

        # Input f_spm needs to be a numpy array
        #if type(f_spm) != np.ndarray:
        if not isinstance(f_spm, np.ndarray):
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
            spm_range = [min(f_spm), max(f_spm)]
            self.__set_sfdu_unpack(spm_range)
        except TypeError as err:
            print('f_spm must be an array of floats or integers: ' +
                '{}'.format(err))
            sys.exit()

        if verbose:
            print('\tAssembling arrays of frequency polynomials from RSR file...')

        # Arrays to contain sets of frequency polynomials
        rfif_lo_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        ddc_lo_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        freq_poly1_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        freq_poly2_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        freq_poly3_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        time_stamp_array = np.zeros(self.__end_sfdu - self.__start_sfdu + 1)
        n_iter = 0
        for i_sfdu in range(self.__start_sfdu, self.__end_sfdu + 1):
            sfdu = self.__rsr_struct[i_sfdu * self.__sfdu_len:
                i_sfdu * self.__sfdu_len + self.__sfdu_len]

            # If end of file reached
            if len(sfdu) == 0:
                break

            # Unpack SFDU into readable format
            s = self.__sfdu_unpack(sfdu)
            s_dict = dict(zip(self.__field_names, s))

            # Beginning of time range for frequency polynomials for this SFDU
            _time_stamp = spm_vals[i_sfdu * self.__n_pts_per_sfdu]

            rfif_lo_array[n_iter] = s_dict['sh_rfif_lo']
            ddc_lo_array[n_iter] = s_dict['sh_ddc_lo']
            freq_poly1_array[n_iter] = s_dict['sh_schan_freq_poly_coef_1']
            freq_poly2_array[n_iter] = s_dict['sh_schan_freq_poly_coef_2']
            freq_poly3_array[n_iter] = s_dict['sh_schan_freq_poly_coef_3']
            time_stamp_array[n_iter] = round(_time_stamp, round_decimal)

            n_iter += 1

        if verbose:
            print('\tEvaluating sky frequency at desired SPM...')

        # Calculate sky frequency for the input time array
        f_sky_pred = np.zeros(len(f_spm))
        for i in range(len(f_spm)):
            # Find correct frequency info for time
            _ind = np.arange(len(time_stamp_array))
            _time_ind = (_ind[f_spm[i] >= time_stamp_array])[-1]

            _msec = f_spm[i] - f_spm[i].astype(int)
            f_sky_pred[i] = ((rfif_lo_array[_time_ind]
                + ddc_lo_array[_time_ind]) * 1.0e6
                - freq_poly1_array[_time_ind]
                - freq_poly2_array[_time_ind] * _msec
                - freq_poly3_array[_time_ind] * (_msec ** 2))

            if verbose and (i < 10):
                print(_msec, f_spm[i], _time_ind, time_stamp_array[_time_ind],
                    freq_poly1_array[_time_ind], freq_poly2_array[_time_ind],
                    freq_poly3_array[_time_ind])

        if verbose:
            for i in range(10):
                print('%24.16f    %30.16f' % (f_spm[i], f_sky_pred[i]))

        # Return f_spm and evaluated predicted sky frequency
        return f_spm, f_sky_pred

    def __set_IQ(self, verbose=False):
        """
        Read full RSR file to find the raw measured I and Q over the
        specified spm_range of the file. Adds attributes for raw SPM values
        over the specified time range and the raw measured I and Q values.

        Arguments
            :verbose (*bool*):
                If True, print steps and intermediate results

        Attributes
            :spm_vals (*np.ndarray*):
                Raw resolution SPM values over specified
                ``spm_range``
            :IQ_m (*np.ndarray*):
                Raw measured complex signal over the specified
                ``spm_range``
        """


        spm_vals = self.spm_vals

        # Record what original SPM values were
        if self.sample_rate_khz == 16:
            self.__spm_16khz = spm_vals


        decimate_16khz_to_1khz = self.__decimate_16khz_to_1khz

        spm_range = [min(spm_vals), max(spm_vals)]
        self.__set_sfdu_unpack(spm_range)

        # Reduce SPM array to match the I and Q arrays to be made
        spm_vals = self.spm_vals[self.__n_pts_per_sfdu * self.__start_sfdu:
            self.__n_pts_per_sfdu * (self.__end_sfdu + 1)]

        # Multiprocessing to retrieve data from RSR file
        results = []
        queues = [Queue() for i in range(self.__cpu_count)]
        n_loops = self.__end_sfdu - self.__start_sfdu + 1
        n_per_core = int(np.floor(n_loops / self.__cpu_count))
        loop_args = [(i * n_per_core, (i + 1) * n_per_core, n_loops,
            queues[i]) for i in range(self.__cpu_count)]
        loop_args[-1] = ((self.__cpu_count - 1) * n_per_core,
            self.__end_sfdu + 1, n_loops, queues[-1])
        jobs = [Process(target=self.__loop, args=(a)) for a in loop_args]
        for j in jobs:
            j.start()
        for q in queues:
            results.append(q.get())
        for j in jobs:
            j.join()
        IQ_m = np.hstack(results)

        # Decimate 16kHz file to 1kHz spacing if specified
        if decimate_16khz_to_1khz & (self.sample_rate_khz == 16):

            if verbose:
                print('\tDecimating to 1kHz sampling...')

            IQ_m = decimate(IQ_m, 4, zero_phase=True)
            IQ_m = decimate(IQ_m, 4, zero_phase=True)

            n_pts = len(IQ_m)
            dt = 1.0 / float(1000)
            spm_vals = spm_vals[0] + dt * np.arange(n_pts)
        elif decimate_16khz_to_1khz & (self.sample_rate_khz == 1) and (
                verbose is True):
            print('\nWARNING (RSRReader.get_IQ): Cannot decimate a 1 kHz file '
                + 'any further. Skipping extra decimation\n')

        self.spm_vals = spm_vals
        self.IQ_m = IQ_m

    def __loop(self, i_start, i_end, n_loops, queue=0):
        """
        Purpose:
            Function to perform loop for multiprocessing

        Arguments:
            :i_start (*int*):
                SFDU number to start indexing at
            :i_end (*int*):
                SFDU number to stop indexing at
            :n_loops (*int*):
                Number of loops that this for loop will go through
                for each processor
            :queue (*object*):
                multiprocessing.Queue instance
        """

        I_array = np.zeros(len(range(i_start, i_end)) * self.__n_pts_per_sfdu)
        Q_array = np.zeros(len(range(i_start, i_end)) * self.__n_pts_per_sfdu)
        i_iter = 0
        for i_sfdu in range(i_start, i_end):
            sfdu = self.__rsr_struct[i_sfdu * self.__sfdu_len:
                i_sfdu * self.__sfdu_len + self.__sfdu_len]

            # If EOF is reached
            if len(sfdu) == 0:
                break

            # Unpack SFDU into readable format
            s = self.__sfdu_unpack(sfdu)
            s_dict = dict(zip(self.__field_names, s))

            I_array[i_iter * self.__n_pts_per_sfdu:
                (i_iter + 1) * self.__n_pts_per_sfdu] = (
                s[-2 * self.__n_pts_per_sfdu + 1::2])
            Q_array[i_iter * self.__n_pts_per_sfdu:
                (i_iter + 1) * self.__n_pts_per_sfdu] = (
                s[-2 * self.__n_pts_per_sfdu::2])

            i_iter += 1

        queue.put(I_array + 1j * Q_array)

    def __set_history(self):
        """
        Purpose:
            Set Python dictionary recording information about the run as the
            history attribute.
        """
        input_var_dict = {'rsr_file': self.rsr_file}
        input_kw_dict = {
            'decimate_16khz_to_1khz': self.__decimate_16khz_to_1khz}

        self.history = write_history_dict(input_var_dict, input_kw_dict,
                __file__)

    def __set_rev_info(self):
        self.rev_info = get_rev_info(self)
"""
Revisions:
"""
