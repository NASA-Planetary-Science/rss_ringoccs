#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
if [ "$#" != 1 ] && [ "$#" != 2 ]; then
    echo "Illegal number of arguments: 2 required. Number given: $#"
    echo "Useage: get_rsr_files.sh list_of_rsr_files.txt /path/to/data"
    exit 1
fi
web="https://atmos.nmsu.edu/pdsd/archive/data/"
input="$1"
while IFS= read -r dirfil
do
    if [ "$#" = 1 ]; then
        localdirfil=$dirfil
    elif [ "$#" = 2 ]; then
        localdirfil="$2$dirfil"
    else
    echo "Illegal number of arguments: 2 required. Number given: $#"
    echo "Useage: get_rsr_files.sh list_of_rsr_files.txt /path/to/data"
    exit 1
        exit 1
    fi
    webdirfil=$web$dirfil
    if [ ! -e $localdirfil ]
    then
        echo Current File $localdirfil does not locally exist. Fetching.
        if curl --output /dev/null --silent --head --fail "$webdirfil"; then
            echo "URL exists: $webdirfil"
            curl --create-dirs -o $localdirfil $webdirfil
        else
            echo "URL does not exist: $webdirfil"
            exit 1
        fi
    else
            echo Current File "$localdirfil" already exists.
    fi
done < "$input"
