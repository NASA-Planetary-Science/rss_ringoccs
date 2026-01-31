#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
    echo "Please use BASH to run this script ($0)" 1>&2
    exit 1
fi

################################################################################
#                                    LICENSE                                   #
################################################################################
#   This file is part of rss_ringoccs.                                         #
#                                                                              #
#   rss_ringoccs is free software: you can redistribute it and/or modify       #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   rss_ringoccs is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with rss_ringoccs.  If not, see <https://www.gnu.org/licenses/>.     #
################################################################################
. _get_rsr_data.sh
InputFile="../tables/rsr_1kHz_files_before_USO_failure.txt"
_get_rsr_data "$InputFile"
