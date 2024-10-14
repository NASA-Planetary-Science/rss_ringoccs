"""
:Purpose:
    Calculate occultation geometry parameters listed within CORSS_8001 v2
    *GEO.TAB and *GEO.LBL files, as well as geometrical quantities needed for
    calibrating RSS ring data.
"""

from .occgeo import Geometry
#from . occgeo_uranus import Uranus_Geometry
from . import calc_occ_geometry
