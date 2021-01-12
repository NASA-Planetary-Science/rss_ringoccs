#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
web="https://atmos.nmsu.edu/pdsd/archive/data/"
# pre-USO failure 1 kHz RSR files
input=../tables/rsr_1kHz_files_before_USO_failure.txt

# Create an error log.
touch get_rsr_1khz_error_log.txt

while IFS= read -r dirfil
do
    localdirfil="../data/$dirfil"
    webdirfil=$web$dirfil
    if [ ! -e $localdirfil ]
    then
        echo Current File $localdirfil does not locally exist. Fetching.
        if curl --output /dev/null --silent --head --fail "$webdirfil"; then
            echo "URL exists: $webdirfil"
            curl --create-dirs -o $localdirfil $webdirfil
        else
            echo "URL does not exist: $webdirfil"
            echo "URL does not exist: $webdirfil" >> get_rsr_1khz_error_log.txt
        fi
    else
            echo Current File "$localdirfil" already exists.
    fi
done < "$input"

# If the error log is empty, remove it.
if [ ! -s get_rsr_1khz_error_log.txt ]
then
    rm -f get_rsr_1khz_error_log.txt
fi
