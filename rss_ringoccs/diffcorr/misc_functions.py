def get_range_request(rngreq):
    """
        Function:
            get_range_request
        Purpose:
            This function takes a variety of input ranges and then
            converts them into the form [start,finish].
        Variables:
            rngreq: The start/end points for diffraction correction.
                    Preferred input is rngreq = [a,b]. Arrays are
                    allowed and the range will be set as:
                        rng = [MIN(array),MAX(array)]
                    Finally, certain strings containing a few of the
                    regions of interests within the rings of Saturn
                    are allowed. Permissible strings are:
                        'maxwell', 'titan', 'huygens', and 'encke'.
                    Strings are neither case nor space sensitive.
                    For other planets use rng = [a,b]
        Outputs:
            rng:    A numpy array of the form [start,finish]
        Dependencies:
            [1] numpy
        Notes:
            [1] Inputs should be in kilometers. Numpy arrays, lists,
                and strings are allowed.
            [2] Strings are neither case nor space sensitive.
            [3] If you enter an illegal string, you will be shown a
                list of allowed strings and prompted for one of them.
            [4] If you enter an illegal value you will be prompted
                for a start and end point.
        References:
            [1] A.F. Cook, F.A. Franklin, F.D. Palluconi,
                Saturn's rings—A survey,
                Icarus, Volume 18, Issue 2, 1973, Pages 317-337,
                https://doi.org/10.1016/0019-1035(73)90214-5
            [2] Pollack, J.B., The Rings of Saturn,
                Space Science Reviews, Volume 18, Issue 1, 1975,
                Pages 3–93, https://doi.org/10.1007/BF00350197
        Examples:
            Try some inputs:
                In [1]: import diffcorr as dc
                In [2]: dc.get_range_request('maxwell')
                Out[2]: array([87410., 87610.])
                In [3]: dc.get_range_request('bob')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: 'titan'
                Out[3]: array([77870., 77930.])
            You don't need to wrap your selection in single quotes.
                In [4]: dc.get_range_request('carl')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: huygens
                Out[4]: array([117650., 117950.])
            You can use caps and double quotes too!
                In [5]: dc.get_range_request('jeff')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: "ENCKE"
                Out[5]: array([132900., 134200.])
            Spacing doesn't matter!
                In [6]: dc.get_range_request('george')
                Illegal range request. Allowed strings are:
                'maxwell', 'titan', 'huygens', 'encke', or 'all'
                Please enter requested range string: m      aX     w el L
                Out[6]: array([87410., 87610.])
                In [7]: dc.get_range_request([123,456])
                Out[7]: array([123., 456.])
                In [8]: dc.get_range_request(6)
                You only provided one range value. I need two.
                Enter starting radius (in Kilometers): 6
                Enter final radius (in Kilometers): 10
                Out[8]: array([ 6., 10.])
        History:
            Translated from IDL: RJM - 2018/05/15 2:03 P.M.
    """
    ti = type(rngreq)
    if (not check_real(rngreq)) and (ti != type([1])) and (ti != type("Hi")):
        print("I don't understand your requested range.")
        start  = input("Enter starting radius (in Kilometers): ")
        start  = float(start)
        finish = input("Enter final radius (in Kilometers): ")
        finish = float(finish)
        rng    = np.array([start,finish])
    elif (ti == type('Hello')):
        reg = rngreq.replace(" ","").lower()
        if (reg in region_dict):
            rng = np.array(region_dict[reg])
        else:
            rngreq = ""
            while(not (rngreq in region_dict)):
                print("Illegal range request.\nAllowed strings are:")
                for key in region_dict: print(key)
                rngreq = input("Please enter requested range string: ")
                rngreq = rngreq.replace(" ","").lower()
                rngreq = rngreq.replace("'","")
                rngreq = rngreq.replace('"',"")
            rng = np.array(region_dict[rngreq])
    elif (np.size(rngreq) < 2):
        print("You only provided one range value. I need two.")
        start  = input("Enter starting radius (in Kilometers): ")
        finish = input("Enter final radius (in Kilometers): ")
        start.replace(" ","").replace("'","").replace('"',"").lower()
        finish.replace(" ","").replace("'","").replace('"',"").lower()
        start  = float(start)
        finish = float(finish)
        rng = np.array([start,finish])
    else:
        rng = np.array([np.float(np.min(rngreq)),np.float(np.max(rngreq))])
    return rng

def get_norm_eq(wtype):
    """
        Function:
            get_norm_eq
        Purpose:
            Compute the Normalized Equivalent Width from a given WTYPE.
        Variables:
            wtype:      The name of the window used for processing.
        Output:
            norm_eq:    The normalized equivalent width of wtype.
        Dependencies:
            [1] numpy
        Notes:
            [1] The normalized equivalent width is pre-computed using
                a window with infinite resolution. That is, the
                intergral is perform on the functions exactly, not
                approximated with Riemann sums. This value is more
                accurate than the method used in the compute_norm_eq
                function, however only a few select windows are
                included. These windows are:
                    [1] Rectangular Window..................'rect'
                    [2] Squared Cosine Window...............'coss'
                    [3] Kaiser-Bessel 2.5 Window............'kb25'
                    [4] Kaiser-Bessel 3.5 Window............'kb35'
                    [5] Modified Kaiser-Bessel 2.5 Window...'kbmd'
            [2] The numerical values that this function returns are
                slightly different than those quoted in the reference
                below. The MTR86 paper only evaluates to two decimals
                whereas we have done double precision to 8 decimals.
            [3] The input must be a string. The function is neither
                case nor space sensitive.
            [4] The normalized equivalent width is a unitless value.
        References:
            [1] Essam A. Marouf, G. Leonard Tyler, Paul A. Rosen,
                Profiling Saturn's rings by radio occultation,
                Icarus, Volume 68, Issue 1, 1986, Pages 120-166,
                https://doi.org/10.1016/0019-1035(86)90078-3
        Examples:
            Test some input strings.
                In [1]: import diffcorr as dc
                In [2]: dc.get_norm_eq('kb25')
                Out[2]: 1.65191895
                In [3]: dc.get_norm_eq('kb35')
                Out[3]: 1.92844639
                In [4]: dc.get_norm_eq('coss')
                Out[4]: 1.5
                In [5]: dc.get_norm_eq('rect')
                Out[5]: 1.0
                In [6]: dc.get_norm_eq('kbmd')
                Out[6]: 1.65994218
            Invalid inputs will prompt you for a new one.
                In [7]: dc.get_norm_eq('bob')
                Invalid window type. Please use one of the following:
                'rect', 'coss', 'kb25', 'kb35', or 'kbmd'
                Please enter a window type: kb25
                Out[7]: 1.65191895
        History:
            Translated from IDL: RJM - 2018/05/15 5:11 P.M.
    """
    if (type(wtype) == type('Hello')):
        wtype = wtype.replace(" ", "").lower()
        if (wtype in func_dict) and (np.size(wtype) == 1):
            norm_eq = func_dict[wtype]["normeq"]
        else:
            print("Invalid window type. Please use one of the following:")
            print("'rect', 'coss', 'kb25', 'kb20', kb35', or 'kbmd'")
            wtype   = input("Please enter a window type: ")
            wtype   = wtype.replace(" ", "").lower()
            wtype   = wtype.replace("'","")
            wtype   = wtype.replace('"',"")
            norm_eq = func_dict[wtype]["normeq"]
    else:
        print("Invalid window type. Please use one of the following:")
        print("'rect', 'coss', 'kb25', 'kb20', kb35', or 'kbmd'")
        wtype   = input("Please enter a window type: ")
        wtype   = wtype.replace(" ", "").lower()
        wtype   = wtype.replace("'","")
        wtype   = wtype.replace('"',"")
        norm_eq = func_dict[wtype]["normeq"]
    return norm_eq

def get_range_actual(rho,rng,w_vals):
    """
        Function:
            get_range_actual
        Purpose:
            Compute the possible allowed range for processing, taking
            into consideration available data (rho) and the requested region.
        Variables:
            RHO:        Radial range of the data.
            RANGE:      Requested start/end points for processing.
            W_KM_VALS:  Window width as a function of ring radius.
        Output:
            START:  The allowed starting point for processing.
            N_USED: The number of points allowed for processing.
        History:
            Translated from IDL: RJM - 2018/05/15 3:19 P.M.
    """
    if (not check_real(rho)):
        sys.exit("Rho must be an array of real numbers")
    if (np.min(rho) < 0.0):
        sys.exit("Rho must be positive")
    if (np.size(rng) != 2):
        sys.exit("Range must have format rng = [a,b]")
    if (not check_pos_real(np.min(rng))):
        sys.exit("Range must be positive")
    if (not check_real(w_vals)):
        sys.exit("w_vals must be real")
    if (np.min(w_vals) < 0.0): 
        sys.exit("w_vals must be positive")
    if (np.min(rng) > np.max(rho)):
        sys.exit("Requested range GREATER than available data.")
    if (np.max(rng) < np.min(rho)):
        sys.exit("Requested range LESS than available data.")
    w_max       = np.max(w_vals)
    rho_min_lim = np.min(rho)+np.ceil(w_max/2.0)
    rho_max_lim = np.max(rho)-np.ceil(w_max/2.0)
    rho_start   = rho[np.min((rho >= np.min(rng)).nonzero())]
    rho_end     = rho[np.max((rho <= np.max(rng)).nonzero())]
    rho_min     = np.max([rho_min_lim,rho_start])
    rho_max     = np.min([rho_max_lim,rho_end])
    start       = int(np.min((rho >= rho_min).nonzero()))
    finish      = int(np.max((rho <= rho_max).nonzero()))
    n_used      = 1 + (finish - start)
    return start, n_used
