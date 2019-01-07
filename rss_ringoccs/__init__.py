"""
    License and Copyright:
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
        with this program. If not, see http://www.gnu.org/licenses/.

        This program is part of the rss_ringoccs repository hosted at
        https://github.com/NASA-Planetary-Science/rss_ringoccs and developed
        with the financial support of NASA's Cassini Mission to Saturn.

    Purpose:
        Provide tools for analysis of ring occultation experiments,
        particularly those pertaining to the `Cassini Radio Science
        Experiment <https://pds-rings.seti.org/cassini/rss/>`_,
        based on methods from [MTR1986]_ and [CRSUG2018]_. This software
        package contains methods for reading and extracting RSR data,
        computing occultation geometry, calibrating RSR data, and
        performing diffraction reconstruction for calibrated data at
        different resolutions. Also included are tools for writing
        and reading the data products output by the software.

    Dependencies:
        #. numpy
        #. spiceypy
        #. scipy
        #. time
        #. sys
        #. subprocess
        #. pdb
        #. matplotlib
        #. mayavi
        #. os
        #. platform

    References:
        .. [MTR1986] Marouf, E., Tyler, L, \& Rosen, P. (1986) 'Profiling Saturn's rings
         by radio occultation', Icarus 68, 120–166.
        .. [CRSUG2018] `Cassini Radio Science User's Guide <https://pds-rings.seti.org/cassini/rss/Cassini%20Radio%20Science%20Users%20Guide%20-%2030%20Sep%202018.pdf>`_.
"""

from . import tools
from . import rsr_reader
from . import occgeo
from . import calibration
from . import diffrec

"""
History:
    Created: Team Cassini - 2018/06/14 2:20 P.M.
    Nov 26 2018 - sflury - doc strings updated to match sphinx formatting
"""
