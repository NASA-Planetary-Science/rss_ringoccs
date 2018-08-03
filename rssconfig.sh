#!/bin/bash
if [ ! "$BASH_VERSION" ] ; then
        echo "Please use BASH to run this script ($0)" 1>&2
        exit 1
fi
osstring=`uname`

if [ osstring = "Darwin" ]; then
	fil="Anaconda3-5.1.0-MacOSX-x86_64.sh"
elif [ osstring = "Linux" ];
	fil="Anaconda3-5.1.0-Linux-x86_64.sh"
else
	echo "Operating System not recognized"
	echo "Only MacOSX and Linux supported"
	echo "Exiting script"
	exit 1
fi

web="https://repo.continuum.io/archive/"
webfil="$web$fil"

if [ ! -e "~/anaconda3" ]; then
	echo "Anaconda3 does not currently exist on your computer."
	echo "Running setup scripts now..."
	# Go to your home directory
	cd ~

	# Use cURL to download Anaconda 5.1.0.
	curl -Ok "$webfil"

	# Run the Anaconda installer.
	bash "filename" -b -p ~/anaconda3

	# Remove the shell script from your home directory.
	rm Anaconda3-5.1.0-MacOSX-x86_64.sh

	# Add Anaconda to your path.
	echo 'export PATH="~/anaconda3/bin:$PATH"' >> ~/.bash_profile

	# Update bash_profile.
	source .bash_profile
else
	echo "You already have Anaconda on your computer."
	echo "Running updates..."
fi

# Update conda
conda update conda

# Install Spiceypy
conda install -c https://conda.anaconda.org/andrewannex spiceypy
