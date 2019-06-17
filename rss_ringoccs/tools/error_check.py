import numpy as np

def check_type(input_var, input_type, input_var_name, f_name):
    # Make sure that input files are strings.
    if not isinstance(input_var, input_type):
        raise TypeError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \r\t%s must have type: %s
                \r\tYour input has type: %s
            """ % (f_name, input_var_name,
                   input_type.__name__, type(input_var).__name__)
        )
    else:
        pass
    
    return

def check_type_and_convert(input_var, input_type, input_var_name, f_name):
    if not isinstance(input_var, input_type):
        try:
            input_var = input_type(input_var)
        except (TypeError, ValueError):
            raise TypeError(
                """
                    \n\r\tError Encountered: rss_ringoccs\n
                    \r\t\t%s\n\n
                    \r\t%s must have type: %s.\n
                    \r\tYour input has type: %s\n
                """
                % (f_name, input_var_name,
                   input_type.__name__, type(input_var).__name__)
            )
    else:
        pass
    
    return input_var

def check_positive(input_var, input_var_name, f_name):
    if (np.min(input_var) <= 0.0):
        raise TypeError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \r\t%s must be be positive.
            """ % (f_name, input_var_name)
        )
    else:
        pass
    
    return

def check_non_negative(input_var, input_var_name, f_name):
    if (np.min(input_var) < 0.0):
        raise TypeError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \r\t%s must be be non-negative.
            """ % (f_name, input_var_name)
        )
    else:
        pass
    
    return

def check_two_pi(input_var, input_var_name, f_name, deg=True):
    if deg:
        max_value = 360.0001
    else:
        max_value = 6.2832

    if (np.max(np.abs(input_var)) > max_value):
        raise TypeError(
            """
                \n\r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \r\t%s must be be less than 2 pi (Radians).
            """ % (f_name, input_var_name)
        )
    else:
        pass
    
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
    try:
        # Check that the requested range is a legal input.
        if (not isinstance(rng, str)) and (not isinstance(rng, list)):
            if (np.size(rng) < 2):
                raise TypeError
            elif (np.min(rng) < 0):
                raise ValueError(
                    """
                        \r\n\tError Encountered:\n
                        \r\t\t%s\n\n
                        \r\tMinimum requested range must be positive\n
                        \r\tYour minimum requested range: %f\n
                    """ % (f_name, np.min(rng))
                )
            else:
                rng = [np.min(rng), np.max(rng)]
        elif isinstance(rng, list):
            # Try converting all elements to floating point numbers.
            if (not all(isinstance(x, float) for x in rng)):
                for i in np.arange(np.size(rng)):
                    rng[i] = float(rng[i])
            else:
                pass

            # Check that there are at least two positive numbers.
            if (np.size(rng) < 2):
                raise TypeError
            elif (np.min(rng) < 0.0):
                raise ValueError(
                    """
                        \r\n\tError Encountered:\n
                        \r\t\t%s\n\n
                        \r\tMinimum requested range must be positive\n
                        \r\tYour minimum requested range: %f\n
                    """ % (f_name, np.min(rng))
                )
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
                \r\n\tError Encountered:\n
                \r\t\t%s\n\n
                \r\trng must be a list or a valid string.\n
                \r\tYour input has type: %s\n
                \r\tSet range=[a,b], where a is the STARTING point\n
                \r\tand b is the ENDING point of reconstruction, or\n
                \r\tuse one of the following valid strings:\n%s
            """ % (f_name, type(rng).__name__, erm)
        )
    
    return rng

def check_psitype(psitype, fname):
    psi_types = ["fresnel", "fresnel3", "fresnel4", "fresnel6",
                 "fresnel8", "full", "ellipse"]

    # Cbeck that psitype is a valid string.
    if not isinstance(psitype, str):
        erm = ""
        for key in psi_types:
            erm = "%s\t\t'%s'\n" % (erm, key)
        raise TypeError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \r\tpsitype must have type: str\n
                \r\tYour input has type: %s\n
                \r\tAllowed strings are:\n%s
            """ % (fname, type(psitype).__name__, erm)
        )
    else:
        # Remove spaces, quotes, and apostrophe's from the psitype.
        psitype = psitype.replace(" ", "").replace("'", "")
        psitype = psitype.replace('"', "").lower()

        # Perform error check, print legal inputs if needed.
        if not (psitype in psi_types):
            erm = ""
            for key in psi_types:
                erm = "%s\t\t'%s'\n" % (erm, key)
            raise TypeError(
                """
                    \r\tError Encountered: rss_ringoccs
                    \r\t\t%s\n
                    \r\tInvalid string for psitype.\n
                    \r\tYour input: %s\n
                    \r\tAllowed strings are:\n%s
                """ % (fname, psitype, erm)
            )
        else:
            pass

    return psitype

def check_lengths(input_var_1, input_var_2, 
                  input_var_name_1, input_var_name_2, function_name):
    if (np.size(input_var_1) != np.size(input_var_2)):
        raise IndexError(
            """
                \r\tError Encountered:
                \r\t\t%s\n
                \r\tThe number of points in %s is not
                \r\tequal to the number of points in %s.
            """ % (function_name, input_var_name_1, input_var_name_2)
        )
    else:
        return

def check_is_real(input_var, input_var_name, fname):
    if not (np.all(np.isreal(input_var))):
        raise TypeError(
            """
                \r\tError Encountered: rss_ringoccs
                \r\t\t%s\n
                \r\t%s must be real valued.
                \r\tYour input has type: %s
            """ % (fname, input_var_name, input_var.dtype)
        )
    else:
        pass

    return
