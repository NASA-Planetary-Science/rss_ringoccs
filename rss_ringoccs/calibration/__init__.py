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

Calibration sub-package of rss_ringoccs

"""

from .calc_freq_offset import calc_freq_offset
from .calibration_class import Calibration
from .freq_offset_fit import FreqOffsetFit
from .norm_diff_class import NormDiff
from .power_normalization import Normalization
