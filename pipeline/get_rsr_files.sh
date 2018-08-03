#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
	echo "Please use BASH to run this script ($0)" 1>&2
	exit 1
fi
web="https://atmos.nmsu.edu/pdsd/archive/data/"
input="$1"
while IFS= read -r dirfil
do
	webdirfil=$web$dirfil
	if [ ! -e $dirfil ]
	then
    		echo Current File $dirfil does not locally exist. Fetching.
		curl --create-dirs -o $dirfil $webdirfil
	else
    		echo Current File "$dirfil" already exists.
	fi
done < "$input"
