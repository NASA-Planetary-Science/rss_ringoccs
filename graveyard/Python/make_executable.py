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
#       Makes a file executable from python on a Unix-like system.             #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Date:   2018/05/16                                                         #
################################################################################
"""

import os

def make_executable(path):
    """
        Function:
            make_executable
        Purpose:
            Makes a file executable from within Python.
        Variables:
            path:
                Path to a file.
        Outputs:
            None.
        Dependencies:
            [1] os
        Notes:
            This routine has only been tested using scripts written
            in Bash, and on standard Unix commands.
    """

    # Current permissions of the file.
    mode = os.stat(path).st_mode

    # Bit mask, in octal, for the file.
    mask = 0o444

    # Flip the bits for the "executable" part of the mask to "1".
    mode |= (mode & mask) >> 2

    # Make the file executable with the unix 'chmod' command.
    os.chmod(path, mode)
