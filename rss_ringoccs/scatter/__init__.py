'''
Purpose:
    Compute the spectrogram of the incoherent signal based on
    properties of the measured signal like sample rate. Spectrogram is
    optionally stacked to improve SNR and output to a .TAB file
    named following the software's output file nomenclature.
'''
from .spectrogram import Scatter
