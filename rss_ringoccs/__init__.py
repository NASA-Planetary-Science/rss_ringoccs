"""
    Package:
        rss_ringoccs
    Purpose:
        Provide tools for analysis of ring occultation experiments,
        particularly those pertaining to the Cassini Radio Science
        Experiment. Mathematical procedures for diffraction
        correction exist, as well as geometrical codes and tools, and
        methods for reading RSR and extracting data. In addition,
        there are several 'utilities' to aid the user in their data
        analysis, as well as mathematical functions and routines for
        modelling both diffraction and diffraction correction based
        on the given geometry.
    Dependencies:
        [01] numpy
        [02] spiceypy
        [03] scipy
        [04] time
        [05] sys
        [06] subprocess
        [07] pdb
        [08] matplotlib
        [09] tkinter
        [10] multiprocessing
        [11] mayavi
        [12] os
        [13] platform
    Subpackages:
    
    References:
        [1] https://github.com/NASA-Planetary-Science/rss_ringoccs/
    History:
        Created: Team Cassini - 2018/06/14 2:20 P.M.
"""

from . import tools
from . import rsr_reader
from . import occgeo
from . import calibration
from . import diffrec