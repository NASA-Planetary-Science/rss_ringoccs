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
#       Executes a shell command from python.                                  #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""

import subprocess

def shell_execute(script):
    """
        Function:
            shell_execute
        Purpose:
            Execute a shell script from within Python.
        Variables:
            script:
                A list containing the path to the shell
                script and variables. For example:
                    script = ['path/to/script','v1',...,'vn']
                Elements must be strings.
        Outputs:
            Process:
                An instance of the Popen class from the subprocess module.
                This contains attributes such as the arguments passed to it,
                and other system details.
        Dependencies:
            [1] subprocess
        Notes:
            This routine has only been tested using scripts written
            in Bash, and on standard Unix commands.
        References:
            [1] https://docs.python.org/3/library/subprocess.html
            [2] https://stackoverflow.com/questions/
                3777301/how-to-call-a-shell-script-from-python-code
    """

    string = ""

    for word in script:
        if not isinstance(word, str):
            raise TypeError("Input must be a list of strings")

        string += word

    return subprocess.Popen([string], shell = True)
