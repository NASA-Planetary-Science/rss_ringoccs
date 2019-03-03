#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
	echo "Please use BASH to run this script ($0)" 1>&2
	exit 1
fi
if [ "$#" != 1 ] && [ "$#" != 2 ]; then
	echo "Illegal number of arguments: 2 required. Number given: $#"
	echo "Useage: get_kernels.sh list_of_kernels.txt /path/to/kernels"
	echo "Example: list_of_kernels.txt containing:"
	echo "naif/CASSINI/kernels/spk/061129R_SCPSE_06292_06307.bsp"
	echo "naif/CASSINI/kernels/spk/061213R_SCPSE_06308_06318.bsp"
	echo "./get_kernels list_of_kernels.txt ../kernels/"
	exit 1
fi
web="https://naif.jpl.nasa.gov/pub/"
input="$1"
while IFS= read -r dirfil
do
	if [ "$#" = 1 ]; then
		localdirfil=$dirfil
	elif [ "$#" = 2 ]; then
		localdirfil="$2$dirfil"
	else
		echo "Illegal number of arguments: Two required"
		echo "Useage: get_kernels.sh list_of_kernels.txt /path/to/kernels"
		echo "Example: list_of_kernels.txt containing:"
		echo "naif/CASSINI/kernels/spk/061129R_SCPSE_06292_06307.bsp"
		echo "naif/CASSINI/kernels/spk/061213R_SCPSE_06308_06318.bsp"
		echo "./get_kernels list_of_kernels.txt ../kernels/"
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
