#!/usr/bin/env python
"""
rsr_reader.py
Purpose: Class to create an instance linked to an RSR file.
Revisions:
      gjs_rsr_reader_v2.py
   Feb 14 2018 - gsteranka - Original version
   Feb 20 2018 - gsteranka - Made edits JWF suggested in her copy of the code.
                             Reading RSR header in __init__, since you always
                             have to do that anyway. Class name changed
                             from ReadRSR to RSRReader. "read_hdr" method now
                             private, since user doesn't need to use
   Feb 20 2018 - gsteranka - DSN output has "DSS-" in front of the number
      gjs_rsr_reader_v3.py
   Feb 28 2018 - gsteranka - Copied from v2. Edited names of methods and
                             edited docstrings. Also changed SPM_vals
                             definition so there's no rounding error
   Mar 02 2018 - gsteranka - Edited so predicted sky frequency and I/Q info
                             are not set as attributes to the objects, but
                             rather returned into variables. Changed so now you
                             can't make two different objects with different
                             attributes set. Also eliminated separate header
                             reading routine and am just doing it in __init__
   Mar 07 2018 - gsteranka - Fixed bug where length of SPM_vals outputted by
                             get_IQ doesn't match length of IQ_m
      rsr_reader.py
   Mar 20 2018 - gsteranka - Copy to official version and remove debug steps
   Mar 21 2018 - gsteranka - Added keyword argument in get_IQ method to
                             decimate 16kHz files to 1kHz sampling if True
"""

import os
import struct
import numpy as np
from scipy.signal import decimate

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
        >>> (f_spm_returned, f_sky_pred) = rsr_inst.get_f_sky_pred(f_spm=f_spm)
        >>>
        >>> # Get raw measured I and Q, and raw SPM, over specified spm_range
        >>> (spm_vals, IQ_m) = rsr_inst.get_IQ(spm_range=spm_range)
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
        'sfdu_length'
    ]
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
        'ph_format_code'
    ]
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
        'sh_reserved2a','sh_reserved2b'
    ]
    __sh_format = 'hh'+'BBh'+'hBB'+'BBcBHccBBbBBBBBHHIBBHHHHH'+22*'d'

    # Data
    __data_field_names = [
        'Data_type',
        'Data_length',
        'Data_QI'
    ]
    __data_header_format = 'HH'

    __field_names = (__sfdu_field_names + __ha_field_names + __ph_field_names
                   + __sh_field_names + __data_field_names)


    def __init__(self,rsr_file, TEST=False):
        """
        Sets full path name of RSR file as an attribute to the instance, and
        reads the header of the RSR file, which sets the RSR header attributes
        spm_vals, doy, year, dsn, band, sample_rate_khz. Also sets private
        attributes for reading the RSR file later in the "get" methods.
        Args:
            rsr_file (str): Full path name of a raw RSR file to read
            TEST (bool): Optional boolean variable which, when set to True,
                prints the header attributes that were set
        """

        self.rsr_file = rsr_file

        # Length of SFDU header, and a function to read the header in the
        # proper format
        struct_hdr_fmt = (self.__endian + self.__sfdu_format + self.__ha_format
                          + self.__ph_format + self.__sh_format
                          + self.__data_header_format)
        struct_hdr_len = struct.calcsize(struct_hdr_fmt)
        struct_unpack_hdr = struct.Struct(struct_hdr_fmt).unpack_from

        # Open first header of SFDU and put into format for unpacking below
        with open(self.rsr_file,'rb') as f:
            sfdu_hdr_raw = f.read(struct_hdr_len)
            f.close()

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

        # Default argment for get_IQ method
        self.__decimate_16khz_to_1khz = False

        if TEST:
            print('First 10 raw SPM:')
            print(self.spm_vals[0:10])
            print('Year, DOY, DSN, band, sample_rate:')
            print(str(self.year) + ', ' + str(self.doy) + ', '
                  + str(self.dsn) + ', ' + str(self.band) + ', '
                  + str(self.sample_rate_khz))


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
                                 self.__n_pts_per_sfdu*end_sfdu
                                 ]

        # Format in which to read rest of RSR file one SFDU at a time
        data_format = (self.__data_header_format
                       + np.int(self.__n_pts_per_sfdu)*'hh')

        # Format to read RSR file in
        rsr_fmt = (self.__endian + self.__sfdu_format + self.__ha_format
                   + self.__ph_format + self.__sh_format + data_format
                   )
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
        Returns:
            f_spm (np.ndarray): Array of SPM values that predicted sky frequency
                was evaluated at.
            f_sky_pred (np.ndarray): Predicted sky frequency, calculated from
                the polynomial coefficients in the RSR file
        """

        spm_vals = self.spm_vals

        # Default 1 second spacing over range of spm_vals
        if f_spm is None:
            f_spm = np.arange(min(spm_vals),max(spm_vals),1.0)

        # Range of SPM values
        spm_range = [min(f_spm),max(f_spm)]

        # SPM is often just a little off from what it's supposed to be, so need
        # to later round to a certain number of decimals. 16 kHz files have
        # more decimals
        if self.sample_rate_khz == 1:
            round_decimal = 4
        else:
            round_decimal = 7

        # Set up unpacking of full RSR
        self.__set_sfdu_unpack(spm_range)

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
                                     i_sfdu*self.__sfdu_len + self.__sfdu_len
                                     ]

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
                              - freq_poly3_array[_time_ind]*(_msec**2)
                              )

        if TEST:
            for i in range(10):
                print('%24.16f    %30.16f' % (f_spm[i],f_sky_pred[i]))

        # Return f_spm and evaluated predicted sky frequency
        return f_spm, f_sky_pred


    def get_IQ(self,spm_range=None, decimate_16khz_to_1khz=None, TEST=False):
        """Read full RSR file to find the raw measured I and Q over the
        specified spm_range of the file. Returns raw SPM values over the
        specified time range and the raw measured I and Q values.
        Args:
            spm_range (list): Optional range of SPM values to read RSR file
                over. 2 element list of the form [minimum, maximum]. Default
                is full range of data
            decimate_16khz_to_1khz (bool): Optional Boolean argument which, if
                set to True, decimates 16kHz files down to 1kHz sampling rate.
                Note that this is a sticky keyword - if you set it to True, it
                will be True for any subsequent calls from the instance until
                you explicitly set it to False. This keyword is linked to the
                private attribute __decimate_16khz_to_1khz
            TEST (bool): If True, print first 10 raw measured I and Q
        Returns:
            spm_vals (np.ndarray): Raw resolution SPM values over specified
                spm_range
            IQ_m (np.ndarray): Raw measured complex signal over the specified
                spm_range
        """

        spm_vals = self.spm_vals

        # See if you defined a new Boolean value. Otherwise, just go with
        # the previous setting (default of False)
        if decimate_16khz_to_1khz is not None:
            self.__decimate_16khz_to_1khz = decimate_16khz_to_1khz
        else:
            decimate_16khz_to_1khz = self.__decimate_16khz_to_1khz

        # Default is whole data set
        if spm_range is None:
            spm_range = [min(spm_vals), max(spm_vals)]

        # Set up unpacking of full RSR
        self.__set_sfdu_unpack(spm_range)

        # Reduce SPM array to match the I and Q arrays to be made
        spm_vals = self.spm_vals[self.__n_pts_per_sfdu*self.__start_sfdu:
                                 self.__n_pts_per_sfdu*(self.__end_sfdu+1)]

        I_array = []
        Q_array = []
        for i_sfdu in range(self.__start_sfdu, self.__end_sfdu+1):
            sfdu = self.__rsr_struct[i_sfdu*self.__sfdu_len:
                                     i_sfdu*self.__sfdu_len + self.__sfdu_len
                                     ]

            # If EOF is reached
            if len(sfdu) == 0:
                break

            # Unpack SFDU into readable format
            s = self.__sfdu_unpack(sfdu)
            s_dict = dict(zip(self.__field_names, s))

            I_array.append(s[-2*self.__n_pts_per_sfdu+1::2])
            Q_array.append(s[-2*self.__n_pts_per_sfdu::2])

        # Raw measured I and Q
        I_m = np.reshape(I_array,-1)
        Q_m = np.reshape(Q_array,-1)
        IQ_m = I_m + 1j*Q_m

        if decimate_16khz_to_1khz & (self.sample_rate_khz == 16):
            IQ_m = decimate(IQ_m, 4, zero_phase=True)
            IQ_m = decimate(IQ_m, 4, zero_phase=True)

            n_pts = len(IQ_m)
            dt = 1.0/float(1000)
            spm_vals = spm_vals[0] + dt*np.arange(n_pts)
        elif decimate_16khz_to_1khz & (self.sample_rate_khz == 1):
            print('WARNING (RSRReader.get_IQ): Cannot decimate a 1 kHz file '
                  +'any further. Skipping extra decimation')

        if TEST:
            print('First 10 SPM, I, and Q:')
            for i in range(10):
                print('%24.16f %15i %15i' %
                      (spm_vals[i], I_m[i], Q_m[i]))

        # Return SPM_vals corresponding to speified spm_range and the raw
        # measured I and Q
        return spm_vals, IQ_m


if __name__ == '__main__':
    pass
