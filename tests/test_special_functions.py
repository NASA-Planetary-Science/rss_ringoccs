import numpy
import scipy
from ..rss_ringoccs.diffrec import special_functions

def p_bessel_J0(x):
    try:
        return scipy.special.jv(0, x)
    except:
        return scipy.special.jv(0, x.astype(numpy.float))

def p_bessel_I0(x):
    try:
        return scipy.special.iv(0, x)
    except:
        return scipy.special.iv(0, x.astype(numpy.float))