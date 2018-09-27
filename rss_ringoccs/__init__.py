"""
Copyright (C) 2018 Team Cassini

This program is free software: you  can redistribute it and/or modify it
under the  terms of the GNU  General Public License as  published by the
Free Software Foundation,  either version 3 of the License,  or (at your
option) any later version.

This  program  is distributed  in  the  hope  that  it will  be  useful,
but WITHOUT  ANY  WARRANTY;  without   even  the  implied  warranty  of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy  of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

This program is part of the rss_ringoccs repository hosted at
https://github.com/NASA-Planetary-Science/rss_ringoccs and developed
with the financial support of NASA's Cassini Mission to Saturn.

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
