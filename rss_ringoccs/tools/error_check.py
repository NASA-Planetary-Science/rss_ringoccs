import numpy

def check_type(input_var, input_type, input_var_name, f_name):
    # Make sure that input files are strings.
    if not isinstance(input_var, input_type):
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\t%s must have type: %s\n\r\tYour input has type: %s\n
            """ % (f_name, input_var_name,
                   input_type.__name__, type(input_var).__name__)
        )
    else:
        return

def check_type_and_convert(input_var, input_type, input_var_name, f_name):
    if not isinstance(input_var, input_type):
        try:
            return input_type(input_var)
        except (TypeError, ValueError):
            raise TypeError(
                """
                \n\r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
                \r\t%s must have type: %s\n\r\tYour input has type: %s\n
                """
                % (f_name, input_var_name,
                   input_type.__name__, type(input_var).__name__)
            )
    else:
        return input_var

def check_positive(input_var, input_var_name, f_name):
    if (numpy.min(input_var) <= 0.0):
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\t%s must be be positive.
            """ % (f_name, input_var_name)
        )
    else:
        return

def check_non_negative(input_var, input_var_name, f_name):
    if (numpy.min(input_var) < 0.0):
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\t%s must be be non-negative.
            """ % (f_name, input_var_name)
        )
    else:
        return

def check_two_pi(input_var, input_var_name, f_name, deg=False):
    if deg:
        max_value = 360.0001
    else:
        max_value = 6.2832

    if (numpy.max(numpy.abs(input_var)) > max_value):
        raise TypeError(
            """
            \n\r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\t%s must be be less than 2 pi (Radians).
            """ % (f_name, input_var_name)
        )
    else:
        return

def check_range_input(rng, f_name):
    region_dict = {
        'all':             [1.0, 400000.0],
        'besselbarnard':   [120210.0, 120330.0],
        'bessel-barnard':  [120210.0, 120330.0],
        'cringripples':    [77690.0, 77760.0],
        'encke':           [132900.0, 134200.0],
        'enckegap':        [132900.0, 134200.0],
        'herschel':        [118100.0, 118380.0],
        'herschelgap':     [118100.0, 118380.0],
        'huygens':         [117650.0, 117950.0],
        'huygensringlet':  [117650.0, 117950.0],
        'janusepimetheus': [96200.0, 96800.0],
        'jeffreys':        [118900.0, 119000.0],
        'jeffreysgap':     [118900.0, 119000.0],
        'kuiper':          [119300.0, 119500.0],
        'kuipergap':       [119300.0, 119500.0],
        'maxwell':         [87410.0, 87610.0],
        'maxwellringlet':  [87410.0, 87610.0],
        'russell':         [118550.0, 118660.0],
        'russellgap':      [118550.0, 118660.0],
        'titan':           [77870.0, 77930.0],
        'titanringlet':    [77870.0, 77930.0]
    }

    # Check that the requested range is a legal input.
    try:
        if (not isinstance(rng, str)) and (not isinstance(rng, list)):
            if (numpy.size(rng) < 2):
                raise TypeError
            elif (numpy.min(rng) < 0):
                raise ValueError
            else:
                rng = [numpy.min(rng), numpy.max(rng)]
        elif isinstance(rng, list):

            # Try converting all elements to floating point numbers.
            if (not all(isinstance(x, float) for x in rng)):
                for i in numpy.arange(numpy.size(rng)):
                    rng[i] = float(rng[i])
            else:
                pass

            # Check that there are at least two positive numbers.
            if (numpy.size(rng) < 2):
                raise TypeError
            elif (numpy.min(rng) < 0.0):
                raise ValueError
            else:
                pass
        elif isinstance(rng, str):
            rng = rng.replace(" ", "").replace("'", "").replace('"', "")
            rng = rng.lower()
            if not (rng in region_dict):
                raise TypeError
            else:
                pass
        else:
            raise TypeError
    except TypeError:
        erm = ""
        for key in region_dict:
            erm = "%s\t\t'%s'\n" % (erm, key)
        raise TypeError(
            """
            \r\n\tError Encountered:\n\r\t\t%s\n
            \r\trng must be a list or a valid string.
            \r\tSet range=[a,b], where a is the STARTING point
            \r\tand b is the ENDING point of reconstruction, or
            \r\tuse one of the following valid strings:\n%s
            """ % (f_name, erm)
        )
    except ValueError:
        raise ValueError(
            """
            \r\n\tError Encountered:\n\r\t\t%s\n
            \r\tMinimum requested range must be positive\n
            \r\tYour minimum requested range: %f\n
            """ % (f_name, numpy.min(rng))
        )
    
    return rng

def check_psitype(psitype, fname):
    psi_types = ["fresnel", "fresnel4", "fresnel6", "fresnel8", "full"]

    # Cbeck that psitype is a valid string.
    try:
        if not isinstance(psitype, str):
            raise TypeError
        else:
            # Remove spaces, quotes, and apostrophe's from the psitype.
            psitype = psitype.replace(" ", "").replace("'", "")
            psitype = psitype.replace('"', "").lower()

            # Perform error check, print legal inputs if needed.
            if not (psitype in psi_types):
                raise TypeError
            else:
                return psitype
    except TypeError:
        erm = ""
        for key in psi_types:
            erm = "%s\t\t'%s'\n" % (erm, key)
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\tpsitype must be a string. Allowed strings are:\n%s
            """ % (fname, erm)
        )

def check_wtype(wtype, fname):
    w_types = ["rect", "coss", "kb20", "kb25",
               "kb35", "kbmd20", "kbmd25", "kbmd35"]

    # Cbeck that wtype is a valid string.
    try:
        if not isinstance(wtype, str):
            raise TypeError
        else:
            # Remove spaces, quotes, and apostrophe's from the wtype.
            wtype = wtype.replace(" ", "").replace("'", "")
            wtype = wtype.replace('"', "").lower()

            # Perform error check, print legal inputs if needed.
            if not (wtype in w_types):
                raise TypeError
            else:
                return wtype
    except TypeError:
        erm = ""
        for key in w_types:
            erm = "%s\t\t'%s'\n" % (erm, key)
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\tpsitype must be a string. Allowed strings are:\n%s
            """ % (fname, erm)
        )

def check_lengths(input_var_1, input_var_2, 
                  input_var_name_1, input_var_name_2, function_name):
    if (numpy.size(input_var_1) != numpy.size(input_var_2)):
        raise IndexError(
            """
            \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\tThe number of points in %s is not
            \r\tequal to the number of points in %s.
            """ % (function_name, input_var_name_1, input_var_name_2)
        )
    else:
        return

def check_is_real(input_var, input_var_name, fname):
    if not (numpy.all(numpy.isreal(input_var))):
        raise TypeError(
            """
            \r\tError Encountered: rss_ringoccs\n\r\t\t%s\n
            \r\t%s must be real valued.
            \r\tYour input has type: %s
            """ % (fname, input_var_name, input_var.dtype)
        )
    else:
        return

def check_model(model, fname):
    model_list = ["deltaimpulse", "fromdata", "leftstraightedge",
                  "rightstraightedge", "ringlet", "gap", "squarewave"]

    # Cbeck that wtype is a valid string.
    try:
        if not isinstance(model, str):
            raise TypeError
        else:
            # Remove spaces, quotes, and apostrophe's from the wtype.
            model = model.replace(" ", "").replace("'", "").replace('"', "")
            model = model.lower()

            # Perform error check, print legal inputs if needed.
            if not (model in model_list):
                raise TypeError
            else:
                return model
    except TypeError:
        erm = ""
        for key in model_list:
            erm = "%s\t\t'%s'\n" % (erm, key)
        raise TypeError(
            """
            \r\tError Encountered:\n\r\t\t%s\n
            \r\tIllegal model requested.\n\r\tYour input: %s
            \r\tAllowed options are:\n%s
            """ % (fname, model, erm)
        )
