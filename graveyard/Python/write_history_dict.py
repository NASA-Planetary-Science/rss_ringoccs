"""
################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or              #
#   modify it under the terms of the GNU General Public License as published   #
#   by the Free Software Foundation, either version 3 of the License, or       #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful             #
#   but WITHOUT ANY WARRANTY# without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
#   Purpose:                                                                   #
#       Create dictionary of processing history for an instance.               #
################################################################################
#   Author: Jolene Fong                                                        #
#   Date:   2018/07/11                                                         #
################################################################################
"""
import time
import os
import platform

def write_history_dict(input_vars, input_kwds, source_file):
    """
        Function:
            write_history_dict
        Purpose:
            This creates a dictionary of processing history for an instance.
        Args:
            input_vars (dict):
                Dictionary of all input variables to the instance.
                For inputs that are instances, use the
                history dictionary from that instance.
            input_kwds (dict):
                Dictionary of all input keywords to the instance.
                For inputs that are instances, use the
                history dictionary from that instance.
            source_file (str):
                Full path to the script used to run the instance.
        Output:
            history (dict):
                Dictionary with keys:
                    "User Name"
                    "Host Name",
                    "Run Date"
                    "Python Version"
                    "Operating System",
                    "Source File"
                    "Input Variables"
                    "Input Keywords"
        Dependencies:
            [1] time
            [2] os
            [3] platform
        History:
            2018/07/11: Jolene Fong
                Original.
            2018/07/27: Jolene Fong
                Added source_files input.
            2024/06/04: Ryan Maguire
                Cleaned up and added to graveyard directory.
    """
    user_name = os.getlogin()
    host_name = os.uname()[1]
    run_date = time.ctime() + ' ' + time.tzname[0]
    python_version = platform.python_version()
    operating_system = os.uname()[0]
    src_dir = source_file.rsplit('/',1)[0] +'/'
    src_file = source_file.split('/')[-1]

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
