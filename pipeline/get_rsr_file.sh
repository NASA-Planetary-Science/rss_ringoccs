#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
web="https://atmos.nmsu.edu/pdsd/archive/data/"
rsr_file_name="$1"
rsr_file_pds_dir="$2"
rsr_file_local_dir="$3"
webdirfil=$web$rsr_file_pds_dir$rsr_file_name
localfil=$rsr_file_local_dir$rsr_file_name
if [ ! -e $localfil ]
then
        echo Current File "$localfil" does not locally exist. Fetching.
        curl --create-dirs -o $localfil $webdirfil
else
        echo Current File "$localfil" already exists.
fi


