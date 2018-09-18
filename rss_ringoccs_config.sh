#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
osstring=`uname`

if [ "$osstring" = "Darwin" ]; then
	fil="Anaconda3-5.1.0-MacOSX-x86_64.sh"
elif [ "$osstring" = "Linux" ]; then
	fil="Anaconda3-5.1.0-Linux-x86_64.sh"
else
	echo "Operating System not recognized"
	echo "Only MacOSX and Linux supported"
	echo "Exiting script"
	exit 1
fi

web="https://repo.continuum.io/archive/"
webfil="$web$fil"

# Go to your home directory
cd ~

# Check to see if anaconda3 exists on your computer.
if [ -d anaconda3 ]; then
	echo "You already have Anaconda on your computer."
	echo "Running updates..."
else
	echo "Anaconda3 does not currently exist on your computer."
	echo "Running setup scripts now..."

	# Use cURL to download Anaconda 5.1.0.
	curl -Ok "$webfil"

	# Run the Anaconda installer.
	bash "$fil" -b -p ~/anaconda3

	# Remove the shell script from your home directory.
	rm -f "$fil"

	# Update bash_profile.
	source .bash_profile
fi

# Update conda
conda update conda

# Install Spiceypy
conda install -c https://conda.anaconda.org/andrewannex spiceypy
