#!/usr/bin/env python
"""

freq_offset.py

Purpose: Calculate frequency offset at a specified time spacing. Returns
         spm_vals and freq_offset_vals at the specified spacing

Revisions:
      gjs_freq_offset.py
   2018 Feb 20 - gsteranka - Original version
   2018 Feb 21 - gsteranka - Added evaluation of f_sky_pred at the "f_SPM"
                             in here
      gjs_freq_offset_v2.py
   2018 Mar 02 - gsterankas - Original version. Edited from past version to
                              be organized as functions instead of in a class
      freq_offset.py
   2018 Mar 20 - gsteranka - Copy to official version and remove debug steps
"""

import numpy as np

def calc_freq_offset(rsr_inst, dt_freq=8.192, spm_range=None, TEST=False):
    """Primary function to run. Calculates frequency offset given a set of
    spm_vals and IQ_m, which can be gotten from RSRReader.get_IQ().

    Example:
        >>> # Get RSR file info
        >>> from rsr_reader import RSRReader
        >>> rsr_inst = RSRReader(rsr_file)
        >>>
        >>> # Calculate frequency offset
        >>> from freq_offset import calc_freq_offset
        >>> f_spm, f_offset = calc_freq_offset(rsr_inst, \
                dt_freq=dt_freq, spm_range=spm_range, TEST=TEST)

    Args:
        spm_vals (np.ndarray): Set of raw resolution SPM values that data was
            recorded at
        IQ_m (np.ndarray): Raw measured complex signal corresponding to
            spm_vals
        dt_freq (float): Optional argument for time spacing at which to extract
            frequency offset. Determines length of FFT to do so. Typically gives
            a power of 2 when it's multiplied by the sample rate in Hz (so for
            1kHz files, typical values would be 0.128 sec, 0.256 sec, 0.512 sec,
            1.024 sec, etc.)
        spm_range (list): Optional 2-element list of range of spm_vals over
            which to extract frequency offset. Default is full range of spm_vals
        TEST (bool): Optional testing argument that, if True, prints requested
            and actual start and end times

    Returns:
        f_spm (np.ndarray): Set of SPM values corresponding to the extracted
            frequency offset at time spacing dt_freq
        f_offset (np.ndarray): Extracted frequency offset at time spacing
            dt_freq
        """

    spm_vals = rsr_inst.spm_vals
    IQ_m = rsr_inst.IQ_m

    if spm_range is None:
        spm_range = [spm_vals[0], spm_vals[-1]]

    dt_raw = spm_vals[1] - spm_vals[0]
    sample_rate_hz = round(1.0 / dt_raw)
    pts_per_fft = int(dt_freq*sample_rate_hz)

    # Range surrounding peak frequency that continuous FFT searches
    #     to refine the peak
    cont_fft_freq_range = 5.0

    # start and end of where to get frequency offset
    start_ind = int(min(np.argwhere(abs(spm_vals - spm_range[0])
        <= dt_raw)))
    end_ind = int(max(np.argwhere(abs(spm_vals - spm_range[1]) <= dt_raw)))

    if TEST:
        print('Requested start and end times:')
        print('%24.16f %24.16f\n' % (spm_vals[start_ind], spm_vals[end_ind]))

    # Eliminate remainder points beyond integer number of pts per
    #     FFT. Go a bit over demanded time range if possible
    remainder_pts = len(spm_vals[start_ind:end_ind]) % pts_per_fft
    if remainder_pts != 0:
        end_ind -= remainder_pts
        end_ind += pts_per_fft
        if end_ind > (len(spm_vals)-1):
            end_ind -= pts_per_fft

    if TEST:
        print('Actual start and end times:')
        print('%24.16f %24.16f\n' % (spm_vals[start_ind], spm_vals[end_ind]))

    # Limit to new start/end times
    spm_vals = spm_vals[start_ind:end_ind]
    IQ_m = IQ_m[start_ind:end_ind]

    # Arrays that will contain final SPM and frequency offset
    n_loops = int(len(spm_vals) / pts_per_fft)
    f_spm = np.zeros(n_loops)
    f_offset = np.zeros(n_loops)

    print('SPM, frequency offset')

    for i in range(n_loops):
        _spm = spm_vals[i*pts_per_fft:(i+1)*pts_per_fft]
        _IQ_m = IQ_m[i*pts_per_fft:(i+1)*pts_per_fft]

        _spm_avg = np.mean(_spm)
        _freq = __find_peak_freq(_spm, _IQ_m, cont_fft_freq_range)

        print('%24.16f %24.16f' % (_spm_avg, _freq))

        f_spm[i] = _spm_avg
        f_offset[i] = _freq

    return f_spm, f_offset


def __find_peak_freq(spm,IQ_m, cont_fft_freq_range):
    """Routine called by calc_freq_offset to find peak frequency. Makes 3
    passes over power spectrum, each at successively finer spacing, to find
    a more accurate frequency offset.

    Args:
        spm (np.ndarray): SPM values for data points over which you're getting
            power spectrum
        IQ_m (np.ndarray): Raw measured complex signal to get the power
            spectrum from
        cont_fft_freq_range (float): Range surrouding peak frequency to look
            for a more refined peak frequency. Set to 5 at top of primary
            function calc_freq_offset

    Returns:
        freq_max_pass3 (float): Frequency of maximum power in the power
            spectrum after the 3rd pass of refining the peak frequency.
            This value is the extracted frequency offset
            """

    # Number of pts per FFT in each pass to refine peak frequency
    cont_fft_num_pass1 = 1000
    cont_fft_num_pass2 = 2000
    cont_fft_num_pass3 = 20

    # Redefine needed variables from __init__
    pts_per_fft = len(spm)
    dt_raw = spm[1] - spm[0]

    # Hamming window
    weight = 0.5*(1.0 - np.cos(2.0*np.pi*np.arange(pts_per_fft)/
        (pts_per_fft-1)))
    IQ_m_weight = IQ_m * weight

    # Take FFT and get the power spectrum
    IQ_m_fft = np.fft.fft(IQ_m_weight)
    power = (np.absolute(IQ_m_fft)**2) / (pts_per_fft**2)

    # Array of frequency values corresponding to "power" array
    freq_half_length = int(pts_per_fft/2)
    freq_min = 1.0 / (pts_per_fft*dt_raw)
    freq1 = np.arange(freq_half_length) * freq_min
    freq2 = (np.arange(freq_half_length) - freq_half_length) * freq_min
    freq = np.concatenate((freq1, freq2))

    # Index of max power, and corresponding frequency
    ind_max = np.argmax(power)
    freq_max = freq[ind_max]

    # Frequency centered at the maximum frequency for first pass of
    # refining max frequency.
    freq_c1 = (freq_max - cont_fft_freq_range +
        np.arange(cont_fft_num_pass1)*
        (2.0*cont_fft_freq_range/cont_fft_num_pass1))
    freq_max_pass1 = __refine_peak_freq(IQ_m_weight, freq_c1, dt_raw)

    # Refine max frequency with a second pass
    freq_c2 = (freq_max_pass1 - 1.0 +
        np.arange(cont_fft_num_pass2)*(2.0/cont_fft_num_pass2))
    freq_max_pass2 = __refine_peak_freq(IQ_m_weight, freq_c2, dt_raw)

    # Refine max frequency with a third pass
    freq_step = 1.0 / cont_fft_num_pass2
    freq_c3 = (freq_max_pass2 - freq_step +
        np.arange(cont_fft_num_pass3)*
        (2.0*freq_step/cont_fft_num_pass3))
    freq_max_pass3 = __refine_peak_freq(IQ_m_weight, freq_c3, dt_raw)

    return freq_max_pass3


def __refine_peak_freq(IQ_m_weight, freq_c, dt):
    """Called by __find_peak_freq to refine peak frequency.

    Args:
        IQ_m_weight (np.ndarray): Weighted raw measured complex signal
        freq_c (np.ndarray): Set of frequency values centered at the peak
            frequency from the last pass
        dt (float): Time spacing of the  IQ_m_weight entries

    Returns:
        freq_max (float): Frequency of maximum power"""

    n = len(IQ_m_weight)
    pow_c = np.zeros(len(freq_c))
    for i in range(len(freq_c)):
        theta = 2.0 * np.pi * freq_c[i] * dt
        zexk = np.exp(-1j*theta*(np.arange(n) + 1))
        sum_zexk = np.sum(IQ_m_weight*zexk) + 1.0
        pow_c[i] = (np.absolute(sum_zexk)**2)/(n**2)

    ind_max = np.argmax(pow_c)
    freq_max = freq_c[ind_max]

    return freq_max


if __name__ == '__main__':
    pass
