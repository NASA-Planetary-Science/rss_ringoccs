#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
	echo "Please use BASH to run this script ($0)" 1>&2
	exit 1
fi
web="https://naif.jpl.nasa.gov/pub/"
input=../tables/complete_list_of_kernels.txt
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
