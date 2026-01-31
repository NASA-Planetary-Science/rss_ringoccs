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
#   Purpose:                                                                   #
#       Short bash script for downloading the RSR data from the PDS.           #
################################################################################
#   Author: Ryan Maguire                                                       #
#   Data:   2025/04/02                                                         #
################################################################################

# Downloads the data listed in a given file. Usage:
#   get_rsr_data /path/to/table.txt
get_rsr_data() {

    website="https://atmos.nmsu.edu/pdsd/archive/data/"

    # pre-USO failure 1 kHz RSR files
    InputFile="$1"
    ErrorLog="./get_rsr_data_error_log.txt"

    # Create an error log.
    touch "$ErrorLog"

    while IFS= read -r directoryFile
    do
        localFile="../data/$directoryFile"
        webFile="$website$directoryFile"

        # If this data has already been downloaded we may skip it.
        if [ ! -e "$localFile" ]; then

            # Check if the URL is valid before trying to download from it.
            validURL=""

            if curl --output /dev/null --silent --head --fail $webFile; then
                curl --create-dirs -o "$localFile" "$webFile"
            else
                echo "URL does not exist: $webFile" >> "$ErrorLog"
            fi
        else
            echo "$localFile already exists."
        fi
    done < "$InputFile"

    # If the error log is empty, remove it.
    if [ ! -s "$ErrorLog" ]; then
        rm -f "$ErrorLog"
    fi
}