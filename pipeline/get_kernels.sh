#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
    echo "Please use BASH to run this script ($0)" 1>&2
    exit 1
fi

web="https://naif.jpl.nasa.gov/pub/"
input="../tables/complete_list_of_kernels.txt"

# Create an error log.
touch get_kernels_error_log.txt

while IFS= read -r dirfil
do
    localdirfil="../kernels/$dirfil"
    webdirfil=$web$dirfil
    if [ ! -e $localdirfil ]; then

        echo "Current File $localdirfil. Fetching..."

        if curl --output /dev/null --silent --head --fail "$webdirfil"; then
            curl --create-dirs -o $localdirfil $webdirfil
        else
            echo "URL does not exist: $webdirfil" >> get_kernels_error_log.txt
        fi
    else
        echo Current File "$localdirfil" already exists.
    fi
done < "$input"

# If the error log is empty, remove it.
if [ ! -s get_kernels_error_log.txt ]; then
    rm -f get_kernels_error_log.txt
fi
