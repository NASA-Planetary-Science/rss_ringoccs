"""
write_history_dict.py

Purpose: Create dictionary of processing history for an instance.

Revisions:
    2018 Jul 11 - jfong - original
"""

import sys
import time
import os
import platform


def write_history_dict(input_vars, input_kwds):
    """
    This creates a dictionary of processing history for an instance.

    Args:
        input_vars (dict):  Dictionary of all input variables to the instance.
                            For inputs that are instances, use the history
                            dictionary from that instance.
        input_kwds (dict):  Dictionary of all input keywords to the instance.
                            For inputs that are instances, use the history
                            dictionary from that instance.

    Output:
        history (dict):     Dictionary with keys:
                            "User Name"
                            "Host Name",
                            "Run Date"
                            "Python Version"
                            "Operating System",
                            "Source File"
                            "Input Variables"
                            "Input Keywords"
    
    Dependencies:
        [1] sys
        [2] time
        [3] os
        [4] platform
    """

    user_name = os.getlogin()
    host_name = os.uname()[1]
    run_date = time.ctime() + ' ' + time.tzname[0]
    python_version = platform.python_version()
    operating_system = os.uname()[0]
    src_dir = __file__.rsplit('/',1)[0] +'/'
    src_file = __file__.split('/')[-1]

    history = {
            "User Name": user_name,
            "Host Name": host_name,
            "Run Date": run_date,
            "Python Version": python_version,
            "Operating System": operating_system,
            "Source Directory": src_dir,
            "Source File": src_file,
            "Input Variables": input_vars,
            "Input Keywords": input_kwds
            }
    return history

