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

echo "Running rss_ringoccs_config Script:"

web="https://repo.continuum.io/archive/"
webfil="$web$fil"

# Go to your home directory
MY_DIRECTORY = $(pwd)
cd ~

# Check to see if anaconda3 exists on your computer.
if [ -d anaconda3 ]; then
	echo -e ' \t ' "You already have Anaconda on your computer."
	echo -e ' \t ' "Running updates..."
else
	echo -e ' \t ' "Anaconda3 does not currently exist on your computer."
	echo -e ' \t ' "Running setup scripts now..."

	# Use cURL to download Anaconda 5.1.0.
	curl -Ok "$webfil"

	# Run the Anaconda installer.
	bash "$fil" -b -p ~/anaconda3

	# Remove the shell script from your home directory.
	rm -f "$fil"
fi

# Check that anaconda3 is included in your path
echo -e ' \t ' "Checking your PATH variable..."
if [[ ":$PATH:" == *"~/anaconda3/bin"* ]]; then
  echo -e ' \t ' "PATH variable is set"
else
	echo -e ' \t ' "Adding ~/anaconda3/bin to your PATH"
	echo 'export "PATH=~/anaconda3/bin:$PATH"' >>~/.bash_profile
fi

# Source .bash_profile to update you path.
echo -e ' \t ' "Sourcing .bash_profile"
echo -e ' \t ' "Make sure you are using Bash when running rss_ringoccs"
source .bash_profile

# Update conda
conda update conda

# Install Spiceypy
conda install -c https://conda.anaconda.org/andrewannex spiceypy

# Source the .bash_profile again after conda updates.
source .bash_profile

# Change back to working directory
cd "$MY_DIRECTORY"