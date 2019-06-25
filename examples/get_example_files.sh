#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
# get kernel files
echo "Getting example kernel files..."
web="https://naif.jpl.nasa.gov/pub/"
input=./Rev007_list_of_kernels.txt
while IFS= read -r dirfil
do
    localdirfil="../kernels/$dirfil"
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
# get rsr files
echo "Getting example 1 kHz data files..."
web="https://atmos.nmsu.edu/pdsd/archive/data/"
input=./Rev007_list_of_rsr_files.txt
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
            exit 1
        fi
    else
            echo Current File "$localdirfil" already exists.
    fi
done < "$input"
