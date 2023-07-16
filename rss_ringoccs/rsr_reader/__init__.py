"""

Purpose:

    Reads in raw RSR data and meta data from the file header. Stores
    pertinant information and raw data as attributes for future use.
    Provides method for predicting the offset frequency from
    spacecraft telemetry at the time of observation.

"""

from .rsr_reader import RSRReader
from .rsr_header import rsr_header
