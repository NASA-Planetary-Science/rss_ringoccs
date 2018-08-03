#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
	echo "Please use BASH to run this script ($0)" 1>&2
	exit 1
fi
if [ "$#" -ne 1 ]; then
	echo "Illegal number of arguments: Only one required"
	echo "Useage: get_kernels.sh list_of_kernels.txt"
	echo "Example: list_of_kernels.txt containing:"
	echo "naif/CASSINI/kernels/spk/061129R_SCPSE_06292_06307.bsp"
	echo "naif/CASSINI/kernels/spk/061213R_SCPSE_06308_06318.bsp"
	echo "./get_kernels list_of_kernels.txt"
	exit 1
fi
web="https://naif.jpl.nasa.gov/pub/"
input="$1"
while IFS= read -r dirfil
do
	webdirfil=$web$dirfil
	if [ ! -e $dirfil ]
	then
    		echo Current File $dirfil does not locally exist. Fetching.
		if curl --output /dev/null --silent --head --fail "$url"; then
 			echo "URL exists: $url"
			curl --create-dirs -o $dirfil $webdirfil
		else
			echo "URL does not exist: $webdirfil"
			exit 1
		fi
	else
    		echo Current File "$dirfil" already exists.
	fi
done < "$input"
