"""
    Part A:
        Begin independent functions. These functions are independent
        of other functions that are within the module and are listed
        first. These are the various window types, variable checks, a
        function dictionary, and some mathematical functions. There
        may be dependence on numpy and other imported modules. There
        is a function for executing bash/shell commands as well.
"""

import numpy as np
import pandas as pd
import spiceypy as spice
from scipy import interpolate
import sys
import subprocess
import time

class __valid_types:
    """
        Class:
            __valid_types
        Purpose:
            Contains listls of legal input types allowed in the
            various functions defined in the diffcorrpy module.
        Variables:
            There are no variables to this class.
        Outputs:
            rtypes: A list of allowed variable types that are 
                    "Real Valued" in nature. Int, Float, and numpy
                    arrays are included.
            ctypes: A list of allowed variable types that are
                    "Complex Valued" in nature. Complex and complex
                    numpy arrays are included.
            atypes: For "All types." This is the concatenation of 
                    rtypes and ctypes.
        Dependencies:
            [1] numpy
        Notes:
            Lists are NOT included in any of the output lists. Some
            functions in diffcorrpy allow lists as inputs, such as
            function pertaining to a requested range, but numpy
            arrays are also allowed. When in doubt, use numpy arrays.
        History:
            Made into a class: RJM - 2018/05/21 7:52 A.M.
        """
    def __init__(self):
        self.rtypes = None
        self.ctypes = None
        self.atypes = None
    
        self.rtypes = [type(0),
            type(0.),
            type(np.intc(0)),
            type(np.intp(0)),
            type(np.int8(0)),
            type(np.int16(0)),
            type(np.int32(0)),
            type(np.int64(0)),
            type(np.array(0)),
            type(np.float16(0.0)),
            type(np.float32(0.0)),
            type(np.float64(0.0)),
            type(np.arange(10))]

        self.ctypes = [type(np.complex(0)),
            type(np.complex64(0)),
            type(np.complex128(0)),
            type(np.array(range(10)))]

        self.atypes = self.rtypes+self.ctypes

vtypes = __valid_types()

def shell_execute(script):
    """
        Function:
            shell_execute
        Purpose:
            Execute a shell script from within Python.
        Variables:
            script:     A list containing the path to the shell
                        script and variables. For example:
                            script = ['path/to/script','v1',...,'vn']
                        Elements must be strings.
        Outputs:
            Process:    An instance of the Popen class from the
                        subprocess module. This contains attributes
                        such as the arguments passed to it, and other
                        system details.
        Dependencies:
            [1] subprocess
            [2] sys
        Notes:
            This routine has only been tested using scripts written
            in Bash, and on standard Unix commands.
        References:
            [1] https://docs.python.org/3/library/subprocess.html
            [2] https://stackoverflow.com/questions/
                3777301/how-to-call-a-shell-script-from-python-code
        Examples:
            Suppose we have the following shell script test.sh:
                #!/bin/bash
                printf "Hello, World! My name is %s!" "$1"
            Run this shell script inside of Python:
                In [1]: import diffcorr as dc
                In [2]: dc.shell_execute(['./test.sh','Bob'])
                        Hello World! My name is Bob!
            We can also execute simple Unix commands.
                In [1]: import diffcorr as dc
                In [2]: a = dc.shell_execute(["echo","Bob"])
                        Bob
                In [3]: a.args
                Out[3]: [' echo Bob']
        History:
            Created: RJM - 2018/05/16 5:49 P.M.
    """
    string = ""
    for x in script:
        if (type(x) != type("Hi")):
            sys.exit("Input must be a list of strings")
        else:
            string = "%s %s" % (string,x)
    Process=subprocess.Popen([string],shell=True)
    return Process

def date_string():
    """
        Function:
            date_string
        Purpose:
            Create string of the form "yyyy_mm_dd_HH_MM_SS_"
        Variables:  There are no variables to this function.
        Outputs:
            date:   Current date "year/month/day/hour/minute/second"
        Dependencies:
            [1] time
        Notes:
            The end of the string has an underscore "_"
        Examples:
            Get the current date.
                In [1]: import diffcorr as dc
                In [2]: dc.date_string()
                Out[2]: '2018_05_27_11_28_22_'
        History:
            Created: RJM - 2018/05/27 11:27 A.M.
    """
    strings = time.strftime("%Y,%m,%d,%H,%M,%S")
    t = strings.split(',')
    date=""
    for x in t:
        date = "%s%s_" % (date,x)
    return date

def check_boole(boo):
    """
        Function:
            check_boole
        Purpose:
            Check if an input is Boolean.
        Variables:
            boo:    Any valid input (Int, float, string, list, etc).
        Outputs:
            bout:   A Boolean (True or False), depending on whether
                    or not boo is Boolean.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the rtypes list, this function
            returns False. Even if the input is a list like [0] or
            [1], which looks Boolean, the function will return False.
            To see the entire rtypes list, do the following:
                In [1]: import diffcorr as dc
                In [2]: dc.vtypes.rtypes
        References:
            [1] https://en.wikipedia.org/wiki/Boolean_data_type
        Examples:
            Check several inputs and see if they're Boolean.
                In [1]: import diffcorr as dc
                In [2]: dc.check_boole(1)
                Out[2]: True
                In [3]: dc.check_boole(2)
                Out[3]: False
                In [4]: dc.check_boole(0)
                Out[4]: True
                In [5]: dc.check_boole(-1)
                Out[5]: False
                In [6]: dc.check_boole(True)
                Out[6]: True
                In [7]: dc.check_boole(False)
                Out[7]: True
                In [8]: dc.check_boole(1.0)
                Out[8]: True
                In [9]: dc.check_boole(1.5)
                Out[9]: False
                In [10]: dc.check_boole(1+1j)
                Out[10]: False
                In [11]: dc.check_boole([1])
                Out[11]: False
                In [12]: dc.check_boole([0,1])
                Out[12]: False
        History:
            Translated from IDL: RJM - 2018/05/14 1:14 P.M.
            Allow True and False inputs: RJM - 2018/05/16 1:15 P.M.
    """
    tb = type(boo)
    if not (tb in vtypes.rtypes) and (tb != type(True)):
        bout = False
    elif (np.size(boo) != 1):
        bout = False
    elif (boo != 0) and (boo != 1) and (boo != True) and (boo != False):
        bout = False
    else:
        bout = True
    return bout

def check_ternary(ter):
    """
        Function:
            check_ternary
        Purpose:
            Check if an input is ternary.
        Variables:
            ter:    Any valid input (Int, float, string, list, etc).
        Outputs:
            tout:   A Boolean (True or False), depending on whether
                    or not boo is Ternary.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the rtypes list, this function
            returns False. Even if the input is a list like [0] or
            [1], which looks Boolean, the function will return False.
            To see the entire rtypes list, do the following:
                In [1]: import diffcorr as dc
                In [2]: dc.vtypes.rtypes
        References:
            [1] https://en.wikipedia.org/wiki/Three-valued_logic
        Examples:
            Check if several types are Ternary:
                In [1]: import diffcorr as dc
                In [2]: dc.check_ternary(1)
                Out[2]: True
                In [3]: dc.check_ternary(0)
                Out[3]: True
                In [4]: dc.check_ternary(2)
                Out[4]: True
                In [5]: dc.check_ternary(3)
                Out[5]: False
                In [6]: dc.check_ternary(2.5)
                Out[6]: False
                In [7]: dc.check_ternary(-1)
                Out[7]: False
                In [8]: dc.check_ternary([0])
                Out[8]: False
                In [9]: dc.check_ternary([0,1])
                Out[9]: False
                In [10]: dc.check_ternary(1+1j)
                Out[10]: False
        History:
            Translated from IDL: RJM - 2018/05/14 1:36 P.M.
    """
    tt = type(ter)
    if not (tt in vtypes.rtypes):
        tout = False
    elif (np.size(ter) != 1):
        tout = False
    elif (ter != 0) and (ter != 1) and (ter != 2):
        tout = False
    else:
        tout = True
    return tout

def check_pos_real(pos):
    """
        Function:
            check_pos_real
        Purpose:
            Check if an input is a SINGLE positive number.
        Variables:
            pos:    Any valid input (Int, float, string, list, etc).
        Outputs:
            pout:   A Boolean (True or False), depending on whether
                    or not pos is a SINGLE positive real number.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the rtypes list, this function
            returns False. Even if the input is a list like [0] or
            [1], which looks positive and real, the function will
            return False. To see the rtypes list, do the following:
                In [1]: import diffcorr as dc
                In [2]: dc.vtypes.rtypes
        References:
            [1] John Stillwell, The Real Numbers:
                An Introduction to Set Theory and Analysis,
                Springer International Publishing, 2013
            [2] https://en.wikipedia.org/wiki/Real_number
        Examples:
            Check if several types are positive real numbers.
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: dc.check_pos_real(0)
                Out[3]: False
                In [4]: dc.check_pos_real(1)
                Out[4]: True
                In [5]: dc.check_pos_real(1.5)
                Out[5]: True
                In [6]: dc.check_pos_real(-2)
                Out[6]: False
                In [7]: dc.check_pos_real(np.pi)
                Out[7]: True
                In [8]: dc.check_pos_real([0])
                Out[8]: False
                In [9]: dc.check_pos_real([2])
                Out[9]: False
                In [10]: dc.check_pos_real(np.array([1,2,np.pi]))
                Out[10]: False
            In the last example, note the np.array([1,2,np.pi]) is a
            legal array, and an array of real numbers. However, 
            check_pos_real only returns True if the input is a single
            real number. That example is an array of three numbers.
        History:
            Translated from IDL: RJM - 2018/05/14 1:30 P.M.
    """
    tp = type(pos)
    if not (tp in vtypes.rtypes):
        pout = False
    elif (np.size(pos) != 1):
        pout = False
    elif (pos <= 0):
        pout = False
    else:
        pout = True
    return pout

def check_real(rea):
    """
        Function:
            check_real
        Purpose:
            Check if an input is real valued.
        Variables:
            rea:    Any valid input (Int, float, string, list, etc).
        Outputs:
            rout:   A Boolean (True or False), depending on whether
                    or not rea is real real valued.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the rtypes list, this function
            returns False. Even if the input is a list like [0] or
            [1], which looks real valued, the function will return
            False. To see the entire rtypes list, do the following:
                In [1]: import diffcorr as dc
                In [2]: dc.vtypes.rtypes
        References:
            [1] John Stillwell, The Real Numbers:
                An Introduction to Set Theory and Analysis,
                Springer International Publishing, 2013
            [2] https://en.wikipedia.org/wiki/Real_number
        Examples:
            Check if several types are real valued:
                In [1]: import diffcorr as dc
                In [2]: import numpy as np
                In [3]: dc.check_real(0)
                Out[3]: True
                In [4]: dc.check_real(-1)
                Out[4]: True
                In [5]: dc.check_real('bob')
                Out[5]: False
                In [6]: dc.check_real(1+1j)
                Out[6]: False
                In [7]: dc.check_real([1,2])
                Out[7]: False
                In [8]: dc.check_real(np.array([1,2,3+1j]))
                Out[8]: False
                In [9]: dc.check_real(np.array([1,2,3+0j]))
                Out[9]: True
                In [10]: dc.check_real(np.array([1,2,np.pi]))
                Out[10]: True
        History:
            Translated from IDL: RJM - 2018/05/15 11:56 P.M.
    """
    tr = type(rea)
    if not (tr in vtypes.rtypes):
        rout = False
    elif (np.size(rea) == 1):
        if not (np.real(rea) == rea): 
            rout = False
        else:
            rout = True
    else:
        if not (np.real(rea) == rea).all():
            rout = False
        else:
            rout = True
    return rout

def check_complex(com):
    """
        Function:
            check_complex
        Purpose:
            Check if an input is purely complex valued.
        Variables:
            com:    Any valid input (Int, float, string, list, etc).
        Outputs:
            cout:   A Boolean (True or False), depending on whether
                    or not rea is real real valued.
        Dependencies:
            [1] numpy
        Notes:
            If the input is not in the ctypes list, this function
            returns False. Even if the input is a list like [1j] or
            [1+1j], which looks real valued, the function will return
            False. To see the entire ctypes list, do the following:
                In [1]: import diffcorr as dc
                In [2]: dc.vtypes.ctypes
        References:
            [1] Paul Nahin, An Imaginary Tale: The Story of i,
                Princeton Publishing, 1998
            [2] https://en.wikipedia.org/wiki/Complex_number
        Examples:
            Check if several types are complex valued.
            In [1]: import diffcorr as dc
            In [2]: import numpy as np
            In [3]: dc.check_complex(0)
            Out[3]: False
            In [4]: dc.check_complex(0+1j)
            Out[4]: True
            In [5]: dc.check_complex(np.pi+1j)
            Out[5]: True
            In [6]: dc.check_complex(np.array([0,1,1+1j]))
            Out[6]: True
            In [7]: dc.check_complex(np.array([0,1,1]))
            Out[7]: False
            In [8]: dc.check_complex('Bob')
            Out[8]: False
        History:
            Translated from IDL: RJM - 2018/05/15 12:08 P.M.
    """
    tr = type(com)
    if not (tr in vtypes.ctypes):
        cout = False
    elif (np.size(com) == 1):
        if (np.imag(com) == 0):
            cout = False
        else:
            cout = True    
    else:
        if (np.imag(com) == 0).all():
            cout = False
        else:
            cout = True
    return cout

class extract_csv_data(object):
    """
        Class:
            csv_extract_data
        Purpose:
            Read three csv files (Geo, Cal, and Tau) and return
            an instance containing all necessary attributes to run
            diffraction correction. This instance can be fed directly
            into the diffraction_correction class.
        Variables:
            Geo:    A string that contains the location of
                    the requested Geo file.
            Cal:    A string that contains the location of
                    the requested Cal file.
            Tau:    A string that contains the location of
                    the requested Tau file.
        Attributes:
            rho_km_vals:        Ring radius, in kilometers.
            phi_rad_vals:       Ring azimuth angle, in radians.
            p_norm_vals:        Normalized raw power.
            phase_rad_vals:     Difracted phase, in radians.
            b_rad_vals:         Ring opening angle, in radians.
            d_km_vals:          Distance from the spacecraft to the
                                ring intercept point, in kilometers.
            f_sky_hz_vals:      The frequency of the recieved
                                singals, in Hertz.
            rho_dot_kms_vals:   The time derivative of the ring
                                radius, in kilometers per second.
    """

    def __init__(self, geodata, caldata, dlpdata, occ,verbose=True):
        if (type(geodata) != type("Hi!")):
            sys.exit("geo file must be a string.")
        if (type(caldata) != type("Hi!")):
            sys.exit("cal file must be a string.")
        if (type(dlpdata) != type("Hi!")):
            sys.exit("dlp file must be a string.")
        if (type(occ) != type("Hi!")):
            sys.exit("occ must be a string.")
        occ = occ.replace(" ","").lower()
        occ = occ.replace("'","")
        occ = occ.replace('"',"")
        if (occ != 'ingress') and (occ != 'egress'):
            sys.exit("occ must be either ingress or egress.")

        self.rho_km_vals        = None
        self.phi_rad_vals       = None
        self.p_norm_vals        = None
        self.phase_rad_vals     = None
        self.b_rad_vals         = None
        self.d_km_vals          = None
        self.f_sky_hz_vals      = None
        self.rho_dot_kms_vals   = None

        if verbose: print("Extracting Geo Data...")
        geo_dat = self.__get_geo(geodata)
        if verbose: print("Geo Data Complete.")
        if verbose: print("Extracting Cal Data...")
        cal_dat = self.__get_cal(caldata)
        if verbose: print("Cal Data Complete.")
        if verbose: print("Extracting DLP Data...")
        dlp_dat = self.__get_dlp(dlpdata)
        if verbose: print("DLP Data Complete")
        if verbose: print("Retrieving Variables...")
        rho_km_vals    = dlp_dat[...,0]
        phi_deg_vals   = dlp_dat[...,4]
        raw_tau_vals   = dlp_dat[...,5]
        phase_deg_vals = dlp_dat[...,6]
        b_deg_vals     = dlp_dat[...,11]
        geo_rho        = geo_dat[...,3]
        geo_d          = geo_dat[...,7]
        geo_drho       = geo_dat[...,8]
        f_sky_raw_vals = cal_dat[...,1]
        if verbose: print("Computing Variables...")
        phi_rad_vals   = phi_deg_vals*spice.rpd()
        phase_rad_vals = phase_deg_vals*spice.rpd()
        b_rad_vals     = b_deg_vals*spice.rpd()
        raw_mu         = np.sin(np.abs(b_rad_vals))
        p_norm_vals    = np.exp(-raw_tau_vals/raw_mu)

        if (occ == 'ingress'):
            crange = (geo_drho < 0.0).nonzero()
        elif (occ == 'egress'):
            crange = (geo_drho > 0.0).nonzero()
        else:
            sys.exit("Illegal occ input")

        if verbose: print("Interpolating Data...")
        geo_rho          = geo_rho[crange]
        geo_drho         = geo_drho[crange]
        geo_d            = geo_d[crange]
        d_km_interp      = interpolate.interp1d(geo_rho,geo_d,kind='linear')
        d_km_vals        = d_km_interp(rho_km_vals)
        rho_dot_interp   = interpolate.interp1d(geo_rho,geo_drho,kind='linear')
        rho_dot_kms_vals = rho_dot_interp(rho_km_vals)
        n_rho_vals       = np.size(rho_km_vals)
        n_f_vals         = np.size(f_sky_raw_vals)
        frange        = np.arange(n_f_vals)
        xrange        = np.arange(n_rho_vals)*(n_f_vals-1.0)/(n_rho_vals-1.0)
        fsky_interp   = interpolate.interp1d(frange,f_sky_raw_vals,kind='linear')
        f_sky_hz_vals = fsky_interp(xrange)
        
        self.rho_km_vals      = rho_km_vals
        self.phi_rad_vals     = phi_rad_vals
        self.p_norm_vals      = p_norm_vals
        self.phase_rad_vals   = phase_rad_vals
        self.b_rad_vals       = b_rad_vals
        self.d_km_vals        = d_km_vals
        self.f_sky_hz_vals    = f_sky_hz_vals
        self.rho_dot_kms_vals = rho_dot_kms_vals

        if verbose: print("Data Extraction Complete.")

    def __get_geo(self,geodata):
        dfg = np.genfromtxt(geodata, delimiter=',')
        return dfg

    def __get_cal(self,caldata):
        dfc = np.genfromtxt(caldata, delimiter=',')
        return dfc

    def __get_dlp(self,dlpdata):
        dfd = np.genfromtxt(dlpdata, delimiter=',')
        return dfd

class pure_csv_reader(object):
    def __init__(self,dat):
        if (type(dat) != type("Hi!")):
            sys.exit("Text file must be a string.")
        df = pd.read_csv(dat)
        self.rho_km_vals      = np.array(df.rho_km_vals)
        self.phase_rad_vals   = np.array(df.phase_rad_vals)
        self.p_norm_vals      = np.array(df.p_norm_vals)
        self.phi_rad_vals     = np.array(df.phi_rad_vals)
        self.b_rad_vals       = np.array(df.b_rad_vals)
        self.f_sky_hz_vals    = np.array(df.f_sky_hz_vals)
        self.d_km_vals        = np.array(df.d_km_vals)
        self.rho_dot_kms_vals = np.array(df.rho_dot_kms_vals)